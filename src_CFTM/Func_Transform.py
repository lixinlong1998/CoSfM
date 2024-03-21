from osgeo import osr
import numpy as np

'''
Note: the default coordinate system of document should be set as PJCS.
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


def transPointCoord(point_coord, source_wkt, target_wkt):
    '''
    point_coord = list[X,Y,Z]
    point_coord_transformed = tuple(X,Y,Z)
    '''
    try:
        source = osr.SpatialReference()
        source.ImportFromWkt(source_wkt)
        target = osr.SpatialReference()
        target.ImportFromWkt(target_wkt)
        transform = osr.CoordinateTransformation(source, target)
        point_coord_transformed = transform.TransformPoint(point_coord[0], point_coord[1], point_coord[2])
    except Exception:
        return list(point_coord_transformed)
    finally:
        return list(point_coord_transformed)


def transPointCoord_CLCS2GCCS(point_coord_CLCS, CoordinateTransform):
    '''
    CLCS -> GCCS
    point_coord_CLCS = [X,Y,Z,*1]
    point_coord_GCCS = [X,Y,Z]
    '''
    point_coord_CLCS = np.asarray([point_coord_CLCS[0], point_coord_CLCS[1], point_coord_CLCS[2], 1])
    M = CoordinateTransform[0]
    point_coord_GCCS = np.matmul(M, point_coord_CLCS)
    return list(point_coord_GCCS)


def transPointCoord_GCCS2CLCS(point_coord_GCCS, CoordinateTransform):
    '''
    GCCS -> CLCS
    point_coord_GCCS = [X,Y,Z,1]
    point_coord_CLCS = [X,Y,Z]
    '''
    point_coord_GCCS = np.asarray([point_coord_GCCS[0], point_coord_GCCS[1], point_coord_GCCS[2], 1])
    M_inv = CoordinateTransform[5]
    point_coord_CLCS = np.matmul(M_inv, point_coord_GCCS)
    return list(point_coord_CLCS)


def transPointCoord_GCCS2GDCS(point_coord_GCCS, CoordinateAttribute):
    '''
    GCCS -> GDCS
    point_coord_GCCS = [X,Y,Z]
    point_coord_GDCS = [lon,lat,height]
    '''
    source_wkt = CoordinateAttribute[2]
    target_wkt = CoordinateAttribute[3]
    point_coord_GDCS = transPointCoord(point_coord_GCCS, source_wkt, target_wkt)
    return list(point_coord_GDCS)


def transPointCoord_GDCS2GCCS(point_coord_GDCS, CoordinateAttribute):
    '''
    GDCS -> GCCS
    point_coord_GDCS = [lon,lat,height]
    point_coord_GCCS = [X,Y,Z]
    '''
    source_wkt = CoordinateAttribute[3]
    target_wkt = CoordinateAttribute[2]
    point_coord_GCCS = transPointCoord(point_coord_GDCS, source_wkt, target_wkt)
    return list(point_coord_GCCS)


def transPointCoord_GCCS2PJCS(point_coord_GCCS, CoordinateAttribute):
    '''
    GCCS -> PJCS
    point_coord_GCCS = [X,Y,Z]
    point_coord_PJCS = [X,Y,height]
    '''
    source_wkt = CoordinateAttribute[2]
    target_wkt = CoordinateAttribute[1]
    point_coord_PJCS = transPointCoord(point_coord_GCCS, source_wkt, target_wkt)
    return list(point_coord_PJCS)


def transPointCoord_PJCS2GCCS(point_coord_PJCS, CoordinateAttribute):
    '''
    PJCS -> GCCS
    point_coord_PJCS = [X,Y,height]
    point_coord_GCCS = [X,Y,Z]
    '''
    source_wkt = CoordinateAttribute[1]
    target_wkt = CoordinateAttribute[2]
    point_coord_GCCS = transPointCoord(point_coord_PJCS, source_wkt, target_wkt)
    return list(point_coord_GCCS)


def transPointCoord_GDCS2PJCS(point_coord_GDCS, CoordinateAttribute):
    '''
    GDCS -> PJCS
    point_coord_GDCS = [lon,lat,height]
    point_coord_PJCS = [X,Y,height]
    '''
    source_wkt = CoordinateAttribute[3]
    target_wkt = CoordinateAttribute[1]
    point_coord_PJCS = transPointCoord(point_coord_GDCS, source_wkt, target_wkt)
    return list(point_coord_PJCS)


def transPointCoord_PJCS2GDCS(point_coord_PJCS, CoordinateAttribute):
    '''
    PJCS -> GDCS
    point_coord_PJCS = [X,Y,height]
    point_coord_GDCS = [lon,lat,height]
    '''
    source_wkt = CoordinateAttribute[1]
    target_wkt = CoordinateAttribute[3]
    point_coord_GDCS = transPointCoord(point_coord_PJCS, source_wkt, target_wkt)
    return list(point_coord_GDCS)


def transPointCoord_CLCS2PJCS(point_coord_PJCS, CoordinateTransform, CoordinateAttribute):
    '''
    CLCS -> PJCS
    point_coord_CLCS = [X,Y,Z,*1]
    point_coord_PJCS = [X,Y,height]
    '''
    point_coord_GCCS = transPointCoord_CLCS2GCCS(point_coord_PJCS, CoordinateTransform)
    point_coord_PJCS = transPointCoord_GCCS2PJCS(point_coord_GCCS, CoordinateAttribute)
    return point_coord_PJCS


def transPointCoord_PJCS2CLCS(point_coord_PJCS, CoordinateTransform, CoordinateAttribute):
    '''
    PJCS -> CLCS
    point_coord_PJCS = [X,Y,height]
    point_coord_CLCS = [X,Y,Z]
    '''
    point_coord_GCCS = transPointCoord_PJCS2GCCS(point_coord_PJCS, CoordinateAttribute)
    point_coord_CLCS = transPointCoord_GCCS2CLCS(point_coord_GCCS, CoordinateTransform)
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
