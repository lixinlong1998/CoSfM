import Metashape
import numpy as np
import math
import src_Metashape.FuncMs_Transform as MsTrans

'''
write result to file
'''


def MetashapeType2numpyType(array):
    if isinstance(array, Metashape.Vector):
        numpy_array = np.empty([array.size, ], dtype=float)
        for i in range(array.size):
            numpy_array[i] = array[i]
    elif isinstance(array, Metashape.Matrix):
        numpy_array = np.empty([array.size[0], array.size[1]], dtype=float)
        for i in range(array.size[0]):
            for j in range(array.size[1]):
                numpy_array[i, j] = array[i, j]
    else:
        raise Exception("[Script]    parameter should be a Metashape.Matrix or Metashape.Vector!")
    return numpy_array


def getBoundary(points, crs, M):
    '''
    input:
        points
        crs
        M
    output:
        Boundary_CLCS = [Points_Coord_Max[0], Points_Coord_Min[0], Points_Coord_Max[1], Points_Coord_Min[1]]
        Boundary_PJCS = [Points_Coord_Max[2], Points_Coord_Min[2], Points_Coord_Max[3], Points_Coord_Min[3]]
    '''
    # get boundary
    points_Coord = []
    for point_id, point in points:
        # Coordinates
        point_coord_CLCS = point.coord
        point_coord_PJCS = MsTrans.transPointCoord_CLCS2PJCS(point.coord, crs, M)  # from CLCS to PJCS
        points_Coord.append([point_coord_CLCS[0], point_coord_CLCS[1], point_coord_CLCS[2],
                             point_coord_PJCS[0], point_coord_PJCS[1], point_coord_PJCS[2]])
    points_Coord = np.asarray(points_Coord)
    points_Coord_Max = points_Coord.max(axis=0)
    points_Coord_Min = points_Coord.min(axis=0)
    # boundary_X_MAX,boundary_X_MIN,boundary_Y_MAX,boundary_Y_MIN
    # right(east),left(west),upper(north),lower(south)
    Boundary_CLCS = [points_Coord_Max[0], points_Coord_Min[0], points_Coord_Max[1], points_Coord_Min[1]]
    Boundary_PJCS = [points_Coord_Max[3], points_Coord_Min[3], points_Coord_Max[4], points_Coord_Min[4]]
    return Boundary_CLCS, Boundary_PJCS


def getGridId(Point_Coord_PJCS, Boundary_PJCS, girdsize):
    # get boundary
    boundary_X_MAX, boundary_X_MIN, boundary_Y_MAX, boundary_Y_MIN = Boundary_PJCS
    Coord_PJCS_X = Point_Coord_PJCS[0]
    Coord_PJCS_Y = Point_Coord_PJCS[1]
    r = int((boundary_Y_MAX - float(Coord_PJCS_Y)) / girdsize)  # row distance to origin(top-left)
    c = int((float(Coord_PJCS_X) - boundary_X_MIN) / girdsize)  # col distance to origin(top-left)
    return [r, c]


def creatCube(Points, Boundary_PJCS, girdsize):
    '''
    creat gird
    output:
        CubesLayers = {grid_id:[]}
        CubesPoints = {grid_id:Point_IdList}
        CubesCorners = {grid_id:[topLeft, topRight, leftBottom, rightBottom, upper, lower]}
    '''
    # get boundary
    boundary_X_MAX, boundary_X_MIN, boundary_Y_MAX, boundary_Y_MIN = Boundary_PJCS
    # creat empty Points Grid
    rows = math.ceil((boundary_Y_MAX - boundary_Y_MIN) / girdsize)
    cols = math.ceil((boundary_X_MAX - boundary_X_MIN) / girdsize)
    print('[Script]    rows, cols:', rows, cols)
    Cubes = np.mat(np.zeros([rows, cols], dtype=np.uintc, order='C'))

    # generate dictionary corresponds to girds
    CubesLayers = {}
    CubesPoints = {}
    CubesCorners = {}
    for r in range(rows):
        for c in range(cols):
            CubesLayers[(r, c)] = []
            CubesPoints[(r, c)] = []
            CubesCorners[(r, c)] = []

    # append point_id to correspond cube
    for point_id, Point in Points.items():
        [r, c] = getGridId(Point[0][1], Boundary_PJCS, girdsize)
        CubesPoints[(r, c)].append(point_id)

    # get coordinates of vertexes of cube
    for grid_id, Point_IdList in CubesPoints.items():
        # skip cell without points
        if len(Point_IdList) == 0:
            continue
        # get plane coordinates of four conners of grid
        r, c = grid_id
        boundary_X_MAX, boundary_X_MIN, boundary_Y_MAX, boundary_Y_MIN = Boundary_PJCS
        topLeft = [boundary_X_MIN + c * girdsize, boundary_Y_MAX - r * girdsize]
        topRight = [boundary_X_MIN + (c + 1) * girdsize, boundary_Y_MAX - r * girdsize]
        leftBottom = [boundary_X_MIN + c * girdsize, boundary_Y_MAX - (r + 1) * girdsize]
        rightBottom = [boundary_X_MIN + (c + 1) * girdsize, boundary_Y_MAX - (r + 1) * girdsize]
        # get the upper and lower coordinate of cube
        Cube_Z = []
        for point_id in Point_IdList:
            Z, sigmaZ = Points[point_id][0][1][2], Points[point_id][2][2]  # in (m) and (mm) respectively
            Cube_Z.append(Z - 3 * sigmaZ / 1000)
            Cube_Z.append(Z + 3 * sigmaZ / 1000)
        upper = max(Cube_Z)
        lower = min(Cube_Z)
        # convert to vertex
        CubesCorners[grid_id] = [topLeft, topRight, leftBottom, rightBottom, upper, lower]
    print('cubes number:', len(CubesPoints))
    return CubesLayers, CubesPoints, CubesCorners
