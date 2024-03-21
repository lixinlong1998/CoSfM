import numpy as np
import math
from osgeo import ogr, osr

'''
write result to file
'''


def cross_ZComponent(a, b):
    '''
    a,b are vectors
    return the z componect of vector from crossing a with b
    '''
    return (a[0] * b[1] - a[1] * b[0])


def isInTriangle(point, Triangle):
    '''
    This function judge if the point is in the triangle
    parameters:
        point = np.array([X,Y])
        Triangle = np.ndarray([[X,Y],[X,Y],[X,Y]])
        , note that the vertex order determines the normal direction of the triangle
    theory:
        If the given point is in the triangle, it should be on the same side of the opposite edge as the vertex.
        In other words, if the point is in the triangle, Z1,z2,z3 should have the same sign.
        note that zi = 0 means the point is located on the edge.
    '''
    PA = point - Triangle[0]
    PB = point - Triangle[1]
    PC = point - Triangle[2]
    z1 = cross_ZComponent(PA, PB)
    z2 = cross_ZComponent(PB, PC)
    z3 = cross_ZComponent(PC, PA)
    return z1 * z2 > 0 and z2 * z3 > 0


def isInPolygon(point, Polygon):
    '''
    This function judge if the point is inside in the Convex Polygon
    parameters:
        point = np.array([X,Y])
        Polygon = np.ndarray([[X,Y],[X,Y],...,[X,Y]])
        , note that the vertex order determines the normal direction of the Polygon
    theory:
        If the given point is in the Polygon
        , it should be on the same side of each edge vector(respect to vertex order)
        In other words, Z component of AB cross AP should have the same sign with the normal direction
                        Z component of BC cross BP should have the same sign with the normal direction
                        ...
                        Z component of NA cross NP should have the same sign with the normal direction.
        The normal direction equal to Z component of AB cross BC.
        note that zi = 0 means the point is located on the edge.
    '''
    AB = Polygon[1] - Polygon[0]
    BC = Polygon[2] - Polygon[1]
    normalDirection = cross_ZComponent(AB, BC)
    for vertex_id, vertex in enumerate(Polygon):
        if vertex_id != len(Polygon) - 1:
            edge = Polygon[vertex_id + 1] - vertex
            vector = point - vertex
            if normalDirection * cross_ZComponent(edge, vector) < 0:
                return False
        else:
            edge = Polygon[0] - vertex
            vector = point - vertex
            if normalDirection * cross_ZComponent(edge, vector) < 0:
                return False
    return True


def encloseCubePolygon(vertexes):
    '''
    This function recover 2D polygon of cube from the projection of vertexes of the cube.
    theory:
        Find the convex polygon to enclose scatter.
    parameters:
        vertexes = np.ndarray([[X,Y],[X,Y],[X,Y],[X,Y],[X,Y],[X,Y],[X,Y],[X,Y]])
        ,where the vertexes are the 2D projection points in the image plane
        with X-axis direct to the right, while y-axis direct down.
    output:
        polygon = np.2darray([[X,Y],[X,Y],...,[X,Y]])
    '''
    # get boundary
    vertexesMax = np.max(vertexes, axis=0)
    vertexesMin = np.min(vertexes, axis=0)
    boundaryUp = vertexesMin[1]
    boundaryRight = vertexesMax[0]
    boundaryDown = vertexesMax[1]
    boundaryLeft = vertexesMin[0]
    # add vertexes located on boundary
    bufferUp = []
    bufferRight = []
    bufferDown = []
    bufferLeft = []
    restVertexes = []
    for vertex_id, vertex in enumerate(vertexes):
        # a flag to represent whether the vertex is located on boundary
        inboundary = 0
        # if vertex dose locate on boundary, it will be append to boundary buffer
        if vertex[1] == boundaryUp:
            bufferUp.append([vertex[0], vertex[1], vertex_id])
            inboundary = 1
        if vertex[0] == boundaryRight:
            bufferRight.append([vertex[0], vertex[1], vertex_id])
            inboundary = 1
        if vertex[1] == boundaryDown:
            bufferDown.append([vertex[0], vertex[1], vertex_id])
            inboundary = 1
        if vertex[0] == boundaryLeft:
            bufferLeft.append([vertex[0], vertex[1], vertex_id])
            inboundary = 1
        # if vertex dose not locate on boundary, it will be append to restVertexes
        if inboundary == 0:
            restVertexes.append([vertex[0], vertex[1], vertex_id])
    # Since vertex may located on the conner of boundarys, the boundary buffers above may contain the same vertex.
    # But more importantly, each boundary buffer is never empty.
    # convert the lists to array
    bufferUp = np.asarray(bufferUp)
    bufferRight = np.asarray(bufferRight)
    bufferDown = np.asarray(bufferDown)
    bufferLeft = np.asarray(bufferLeft)
    restVertexes = np.asarray(restVertexes)
    # rank each buffer with clockwise order
    bufferUp_Rank = bufferUp[bufferUp[:, 0].argsort()]  # sort x from min to max
    bufferRight_Rank = bufferRight[bufferRight[:, 1].argsort()]  # sort y from min to max
    bufferDown_Rank = bufferDown[bufferDown[:, 0].argsort()[::-1]]  # sort x from max to min
    bufferLeft_Rank = bufferLeft[bufferLeft[:, 1].argsort()[::-1]]  # sort y from max to min

    # check whether triangle exist in each conner and whether contain a vertex in it.
    if not (bufferUp_Rank[0][0] == boundaryLeft or bufferLeft_Rank[-1][1] == boundaryUp):
        # print('top-left triangle dose not exist!')
        Triangle = np.array([[boundaryLeft, boundaryUp],
                             [bufferUp_Rank[0][0], boundaryUp], [boundaryLeft, bufferLeft_Rank[-1][1]]])
        for vertex in restVertexes:
            if isInTriangle(vertex[:2], Triangle):
                bufferUp_Rank = np.vstack([np.array([vertex[0], vertex[1], vertex[2]]), bufferUp_Rank])
                # print('find one vertex in top-left triangle')
    if not (bufferUp_Rank[-1][0] == boundaryRight or bufferRight[0][1] == boundaryUp):
        # print('top-right triangle dose not exist!')
        Triangle = np.array([[boundaryRight, boundaryUp],
                             [boundaryRight, bufferRight[0][1]], [bufferUp_Rank[-1][0], boundaryUp]])
        for vertex in restVertexes:
            if isInTriangle(vertex[:2], Triangle):
                bufferRight_Rank = np.vstack([np.array([vertex[0], vertex[1], vertex[2]]), bufferRight_Rank])
                # print('find one vertex in top-right triangle')
    if not (bufferRight[-1][1] == boundaryDown or bufferDown_Rank[0][0] == boundaryRight):
        # print('right-bottom triangle dose not exist!')
        Triangle = np.array([[boundaryRight, boundaryDown],
                             [bufferDown_Rank[0][0], boundaryDown], [boundaryRight, bufferRight[-1][1]]])
        for vertex in restVertexes:
            if isInTriangle(vertex[:2], Triangle):
                bufferDown_Rank = np.vstack([np.array([vertex[0], vertex[1], vertex[2]]), bufferDown_Rank])
                # print('find one vertex in right-bottom triangle')
    if not (bufferDown_Rank[-1][0] == boundaryLeft or bufferLeft_Rank[0][1] == boundaryDown):
        # print('left-bottom triangle dose not exist!')
        Triangle = np.array([[boundaryLeft, boundaryDown],
                             [boundaryLeft, bufferLeft_Rank[0][1]], [bufferDown_Rank[-1][0], boundaryDown]])
        for vertex in restVertexes:
            if isInTriangle(vertex[:2], Triangle):
                bufferLeft_Rank = np.vstack([np.array([vertex[0], vertex[1], vertex[2]]), bufferLeft_Rank])
                # print('find one vertex in left-bottom triangle')

    # combine these buffer to get the polygon vertex
    buffer = np.vstack([bufferUp_Rank, bufferRight_Rank, bufferDown_Rank, bufferLeft_Rank])
    polygon = []
    for vertex in buffer:
        vertex_coord = (vertex[0], vertex[1])
        if vertex_coord in polygon:
            continue
        else:
            polygon.append(vertex_coord)
    return np.asarray(polygon)


def creatGrid(Boundary_PJCS, girdsize):
    '''
    creat empty gird with its plane coordinates
    Grid = [row,col,topLeft, topRight, leftBottom, rightBottom]
    '''
    # get boundary
    boundary_X_MAX, boundary_X_MIN, boundary_Y_MAX, boundary_Y_MIN = Boundary_PJCS
    # creat empty PointsGrid and ICTPsGrid
    rows = math.ceil((boundary_Y_MAX - boundary_Y_MIN) / girdsize)
    cols = math.ceil((boundary_X_MAX - boundary_X_MIN) / girdsize)
    print('[Script]    rows, cols:', rows, cols)
    # generate points dictionary corresponds to girds
    Grids = {}
    for r in range(rows):
        for c in range(cols):
            topLeft = [boundary_X_MIN + c * girdsize, boundary_Y_MAX - r * girdsize]
            topRight = [boundary_X_MIN + (c + 1) * girdsize, boundary_Y_MAX - r * girdsize]
            leftBottom = [boundary_X_MIN + c * girdsize, boundary_Y_MAX - (r + 1) * girdsize]
            rightBottom = [boundary_X_MIN + (c + 1) * girdsize, boundary_Y_MAX - (r + 1) * girdsize]
            Grids[(r, c)] = [r, c, topLeft, topRight, leftBottom, rightBottom]
    return Grids


def getRegidTransformation(points):
    '''
    假设你有一个形状为(1000, 6)的NumPy数组points，将连接点数据分割为源点（source points）和目标点（target points）。假设前3列为源点，后3列为目标点。
    '''
    # 分割源点和目标点
    source_points = points[:, :3]
    target_points = points[:, 3:]
    # 计算质心
    centroid_source = np.mean(source_points, axis=0)
    centroid_target = np.mean(target_points, axis=0)
    # 计算去中心化的点对
    source_points_centered = source_points - centroid_source
    target_points_centered = target_points - centroid_target

    # 计算协方差矩阵
    covariance_matrix = np.dot(source_points_centered.T, target_points_centered)
    # 使用奇异值分解（SVD）求解旋转矩阵
    U, _, Vt = np.linalg.svd(covariance_matrix)
    R = np.dot(U, Vt)
    # 计算平移向量
    T = centroid_target - np.dot(R, centroid_source)

    # 计算配准后的三维残差矢量
    registered_source_points = np.dot(R, source_points.T).T + T
    residuals = registered_source_points - target_points

    # print("旋转矩阵R:", R)
    # print("平移向量T:", T)
    # print("配准后的连接点对的残差:", residuals)
    return R, T, residuals

def isInPolygonShp(point, Polygon):
    '''
    This function judge if the point is inside in the Polygon.shp
    parameters:
        point = np.array([X,Y])
        Polygon = shapefile.GetLayer()
    '''
    layer = Polygon.GetLayer()
    # 获取.shp文件的空间参考信息
    spatial_ref = layer.GetSpatialRef()

    # 创建待检查的点
    wkbpoint = ogr.Geometry(ogr.wkbPoint)
    wkbpoint.AddPoint(point[0], point[1])  # x和y是待检查点的二维地理参考坐标,此处假设二维点的坐标系与面要素的坐标系一致
    wkbpoint.AssignSpatialReference(spatial_ref)  # 设置点的空间参考

    # 使用GDAL的几何体的`Intersects`方法检查点是否在.shp文件的范围内
    layer.ResetReading()
    for feature in layer:
        if feature.geometry().Intersects(wkbpoint):
            return True
            break
    else:
        return False

