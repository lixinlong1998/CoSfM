import Metashape
import numpy as np

'''
Metashape actually have four coordinate system types, which are only mentioned three in 'Metashape Python API Document v1.5.5'
1. Chunk local coordinate system, CLCS
   , with origin in region center and scale for numerical calculation.
2. Geocentric coordinate system, GCCS
   , with origin in the earth center and three axis fixed on ellipsoid of the earth.
3. Geographic coordinate system(geodetic or projected), GGCS(GDCS or PJCS)
   , such as WGS84(for geodetic) or WGS84 UTM Zone 47N(for projected).
4. Local space euclidean coordinate system, LSECS
   , with origin in given point and three axis direct to East, North, and Up respectively.
5. camera coordinate system, CAMCS
6. image coordinate system, IMACS
'''


def getCoordTransMat(chunk):
    '''
    M transform matix could transform point from 'chunk local coordinate system' to 'geocentric coordinate'
    M_inv = chunk.transform.matrix.inv()
    '''
    M = chunk.transform.matrix
    T = chunk.crs.localframe(M.mulp(chunk.region.center)) * M
    if chunk.transform.scale:
        R = chunk.transform.scale * T.rotation()
    else:
        R = T.rotation()
    return M, T, R


def transPointCoord_CLCS2GCCS(point_coord_CLCS, M):
    '''
    CLCS -> GCCS
    point_coord_CLCS = Vector([X,Y,Z,*1])
    point_coord_GCCS = Vector([X,Y,Z])
    '''
    point_coord_CLCS = Metashape.Vector([point_coord_CLCS[0], point_coord_CLCS[1], point_coord_CLCS[2]])
    point_coord_GCCS = M.mulp(point_coord_CLCS)
    return point_coord_GCCS


def transPointCoord_GCCS2CLCS(point_coord_GCCS, M):
    '''
    GCCS -> CLCS
    point_coord_GCCS = Vector([X,Y,Z])
    point_coord_CLCS = Vector([X,Y,Z])
    '''
    point_coord_CLCS = M.inv().mulp(point_coord_GCCS)
    return point_coord_CLCS


def transPointCoord_GCCS2GDCS(point_coord_GCCS, crs):
    '''
    GCCS -> GDCS
    point_coord_GCCS = Vector([X,Y,Z])
    point_coord_GDCS = Vector([lon,lat,height])
    '''
    point_coord_GDCS = crs.transform(point_coord_GCCS, crs.geoccs, crs.geogcs)
    return point_coord_GDCS


def transPointCoord_GDCS2GCCS(point_coord_GDCS, crs):
    '''
    GDCS -> GCCS
    point_coord_GDCS = [lon,lat,height]
    point_coord_GCCS = [X,Y,Z]
    '''
    point_coord_GCCS = crs.transform(point_coord_GDCS, crs.geogcs, crs.geoccs)
    return point_coord_GCCS


def transPointCoord_GCCS2PJCS(point_coord_GCCS, crs):
    '''
    GCCS -> PJCS
    point_coord_GCCS = Vector([X,Y,Z])
    point_coord_PJCS = Vector([X,Y,height])
    '''
    point_coord_PJCS = crs.project(point_coord_GCCS)
    return point_coord_PJCS


def transPointCoord_PJCS2GCCS(point_coord_PJCS, crs):
    '''
    PJCS -> GCCS
    point_coord_PJCS = [X,Y,height]
    point_coord_GCCS = [X,Y,Z]
    '''
    point_coord_PJCS = Metashape.Vector([point_coord_PJCS[0], point_coord_PJCS[1], point_coord_PJCS[2]])
    point_coord_GCCS = crs.unproject(point_coord_PJCS)
    return point_coord_GCCS


def transPointCoord_GDCS2PJCS(point_coord_GDCS, crs):
    '''
    GDCS -> PJCS
    point_coord_GDCS = [lon,lat,height]
    point_coord_PJCS = [X,Y,height]
    '''
    point_coord_GCCS = transPointCoord_GDCS2GCCS(point_coord_GDCS, crs)
    point_coord_PJCS = crs.project(point_coord_GCCS)
    return point_coord_PJCS


def transPointCoord_PJCS2GDCS(point_coord_PJCS, crs):
    '''
    PJCS -> GDCS
    point_coord_PJCS = [X,Y,height]
    point_coord_GDCS = [lon,lat,height]
    '''
    point_coord_PJCS = Metashape.Vector([point_coord_PJCS[0], point_coord_PJCS[1], point_coord_PJCS[2]])
    point_coord_GCCS = crs.unproject(point_coord_PJCS)
    point_coord_GDCS = transPointCoord_GCCS2GDCS(point_coord_GCCS, crs)
    return point_coord_GDCS


def transPointCoord_CLCS2PJCS(point_coord_CLCS, crs, M):
    '''
    CLCS -> PJCS
    point_coord_CLCS = [X,Y,Z,*1]
    point_coord_PJCS = Vector([X,Y,height])
    '''
    point_coord_CLCS = Metashape.Vector([point_coord_CLCS[0], point_coord_CLCS[1], point_coord_CLCS[2]])
    point_coord_GCCS = transPointCoord_CLCS2GCCS(point_coord_CLCS, M)
    point_coord_PJCS = transPointCoord_GCCS2PJCS(point_coord_GCCS, crs)
    return point_coord_PJCS



def transPointCoord_PJCS2CLCS(point_coord_PJCS, crs, M):
    '''
    PJCS -> CLCS
    point_coord_PJCS = [X,Y,height]
    point_coord_CLCS = [X,Y,Z,*1]
    '''
    point_coord_PJCS = Metashape.Vector([point_coord_PJCS[0], point_coord_PJCS[1], point_coord_PJCS[2]])
    point_coord_GCCS = transPointCoord_PJCS2GCCS(point_coord_PJCS, crs)
    point_coord_CLCS = transPointCoord_GCCS2CLCS(point_coord_GCCS, M)
    return point_coord_CLCS


def transPoint_CLCS2CAMCS(Point_coord_CLCS, Camera):
    '''
    [X, Y, Z, *1] -> [x_nol, y_nol]
    '''
    Point_coord_CLCS = np.asarray([Point_coord_CLCS[0], Point_coord_CLCS[1], Point_coord_CLCS[2], 1])
    cameraTransformMatrix = Camera[0][0]  # 4darray
    cameraTransformMatrix_inv = np.linalg.inv(cameraTransformMatrix)
    point_CAMCS = np.matmul(cameraTransformMatrix_inv, Point_coord_CLCS)
    point_CAMCS_normal = point_CAMCS[:2] / point_CAMCS[2]
    return [point_CAMCS_normal[0], point_CAMCS_normal[1]]


def transPointCov_CLCS2PJCS(point_cov_CLCS, R):
    '''
    point_cov_CLCS = point.cov
    '''
    # Transform point covariance matrix into the 'output coordinate system'
    point_cov_PJCS = R * point_cov_CLCS * R.t()
    return point_cov_PJCS


def transCameraPose_MAT2XYZOPK(camera_transform, crs, M):
    '''
    convert the camera.transform (in CLCS) to [X,Y,Z,O,P,K](in PJCS & LSECS)
    '''
    # X,Y,Z is in the PJCS
    camera_center = camera_transform.translation()
    xyz = crs.project(M.mulp(camera_center))
    # rotaion angle is in the LSECS(camera center based)
    transform = crs.localframe(M.mulp(camera_center)) * M * camera_transform
    opk = Metashape.Utils.mat2opk(transform.rotation() * Metashape.Matrix.Diag((1, -1, -1)))
    return [xyz[0], xyz[1], xyz[2], opk[0], opk[1], opk[2]]


def transCameraPose_MAT2XYZYPR(camera_transform, crs, M):
    '''
    convert the camera.transform (in CLCS) to [X,Y,Z,Y,P,R](in PJCS & LSECS)
    '''
    # X,Y,Z is in the PJCS
    camera_center = camera_transform.translation()
    xyz = crs.project(M.mulp(camera_center))
    # rotaion angle is in the LSECS(camera center based)
    transform = crs.localframe(M.mulp(camera_center)) * M * camera_transform
    ypr = Metashape.Utils.mat2ypr(transform.rotation() * Metashape.Matrix.Diag((1, -1, -1)))
    return [xyz[0], xyz[1], xyz[2], ypr[0], ypr[1], ypr[2]]
