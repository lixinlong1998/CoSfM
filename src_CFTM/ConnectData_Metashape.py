from config import *
import Metashape
import math
import numpy as np
from src_Metashape import FuncMs_Transform as MsTrans
from src_Metashape import FuncMs_Camera as MsCam
from src_Metashape import FuncMs_Point as Poi


'''
Metashape Python API dose not provide public classes, which makes it unusable for multiprocess module
So, this script is aim at converting data(cameras, points, sensors, projections,etc) from Metashape class to normal one.
The need-converted data are as follows:

    X parameters:   Cameras = {identifier:[q, i, j, k, X, Y, Z, sensor_id],...}
                    Points = {point_id:[X(m), Y(m), Z(m), stdX(mm), stdY(mm), stdZ(mm), RE, RU, PA, IC]}
    
    Coefficient:    Sensors = {sensor_id:[f,cx,cy,b1,b2,k1,k2,k3,k4,p1,p2,p3,p4]}
                    Projections = {identifier:[[u,v,size,track_id],...],...}
                    Tracks = [Track[view,...]]
                    CoordinateTransform = [M,LSE,T,R,S,M_inv]
                    CoordinateAttribute = [...]

    
    Indices:        camera_ids
                    point_ids
                    validTracks
                    validTracks_ids
                    pair_ids
                    match_ids
    
    
                    points2D = [[u,v],...] with shape (nobservations,2)
                    - the observation array
                    new_cameras_id = [camera_id,...] with shape (nnew_cameras,) 
                    - the index of new camera in original cameras
                    camerasModel_id = [sensor_id,...] with shape (nnew_cameras,) 
                    - the index of new camera's sensor in original sensors
                    cameraIndices = [new_camera_id,...] with shape (nobservations,)
                    - the new camera index of observation
                    pointIndices = [new_point_id,...] with shape (nobservations,)
                    - the new point index of observation
'''


def MetashapeType2numpyType(array):
    '''
    convert Metashape.Matrix or Metashape.Vector class to numpy ndarrray
    '''
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
        numpy_array = np.asarray(array)
    return numpy_array


def getComponent_Cameras(chunk, Covariance=True):
    Functime = time.perf_counter()
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)
    Cameras = {}
    for camera_id, camera in enumerate(chunk.cameras):
        # skipping Keyframes
        if camera.type == Metashape.Camera.Type.Keyframe:
            print('[Script]    the sensor is keyframe!:', camera.label)
            continue
        elif not camera.transform:
            print('[Script]    camera has no transform!:', camera.label)
            continue  # skipping no transform
        elif not camera.enabled:
            print('[Script]    camera is not enabled!', camera.label)
            continue  # skipping not enabled

        # construct key
        identifier = MsCam.getCameraIdentifier(camera)

        # construct value
        '''
        [transform, reference, covariance, sensor_id, state,path]
            camera_transform = Cameras[identifier][0][0]
            camera_loc = Cameras[identifier][0][1]
            camera_rot = Cameras[identifier][0][2]
            camera_transform_PJCS = Cameras[identifier][0][3]
            camera_loc_Ref = Cameras[identifier][1][0]
            camera_loc_Acc = Cameras[identifier][1][1]
            camera_rot_Ref = Cameras[identifier][1][2]
            camera_rot_Acc = Cameras[identifier][1][3]
        note:point in camera coordinate left multiplied by camera.transform will transform the point to chunk coordinate
        '''
        ## transform = [camera_transform,camera_loc,camera_rot,camera_transform_PJCS]
        camera_transform = MetashapeType2numpyType(camera.transform)  # return matrix with shape 4x4
        camera_loc = MetashapeType2numpyType(camera.transform.translation())  # return ndarray([x,y,z])
        camera_rot = MetashapeType2numpyType(camera.transform.rotation())  # return matrix with shape 3x3
        # camera_qua = Rotation.from_matrix(camera_rot).as_quat()  # return [q,i,j,k]
        camera_transform_PJCS = MetashapeType2numpyType(
            MsCam.transCameraPose_MAT2XYZOPK(camera.transform, crs, M))  # ndarray([X,Y,Z,O,P,K]) in PJCS&LSE
        transform = [camera_transform, camera_loc, camera_rot, camera_transform_PJCS]
        ## reference = [camera_loc_Ref,location_accuracy,camera_rot_Ref,rotation_accuracy]
        camera_loc_Ref = MetashapeType2numpyType(camera.reference.location)  # return ndarray([x,y,z])
        camera_loc_Acc = MetashapeType2numpyType(camera.reference.location_accuracy)  # return ndarray([x,y,z]) or None
        camera_rot_Ref = MetashapeType2numpyType(camera.reference.rotation)  # return ndarray([x,y,z])
        camera_rot_Acc = MetashapeType2numpyType(camera.reference.rotation_accuracy)  # return ndarray([x,y,z]) or None
        reference = [camera_loc_Ref, camera_loc_Acc, camera_rot_Ref, camera_rot_Acc]
        ## covariance
        if Covariance:
            covariance = [MetashapeType2numpyType(camera.location_covariance),
                          MetashapeType2numpyType(camera.rotation_covariance)]
        else:
            covariance = []
        ## sensor_id
        sensor_id = chunk.sensors.index(camera.sensor)
        ## state =[]
        camera_label = camera.label
        camera_loc_Ref = camera.reference.location_enabled
        camera_rot_Ref = camera.reference.rotation_enabled
        camera_selected = camera.selected
        camera_orientation = camera.orientation
        states = [camera_label, camera_loc_Ref, camera_rot_Ref, camera_selected, camera_orientation]
        ## path = camera.photo.path
        Camera = [transform, reference, covariance, sensor_id, states, camera.photo.path, identifier]
        Cameras[camera_id] = Camera
    print('[Script][TimeCost]    DataImporting_Metashape.getComponent_Cameras:', time.perf_counter() - Functime, 'sec')
    return Cameras


def getComponent_Points(chunk, WriteCovariance=True, WriteQuality=True):
    '''
    get points coordinate, covariance, and quality as required
    output:
        Points = [Point,Point,...,Point]
        Point = [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, track_id]
        Coordinates=[[X, Y, Z],[X, Y, Z]]   # in CLCS & PJCS respectively
        Covariance = 3x3 covariance matrix
        StandardDeviation = [stdX, stdY, stdZ]  # in mm
        Quality = [RE, RU, PA, IC]
    '''
    Functime = time.perf_counter()
    points = chunk.point_cloud.points
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)
    if WriteQuality:
        PointsCriterion = Poi.getPointsCriterion(chunk)

    Points = {}
    for point_id, point in enumerate(points):
        # skipping invalid points
        if not point.valid:
            continue

        # [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, track_id]
        Point = [[], [], [], [], 0, 0, 0]

        # Coordinates
        point_coord_CLCS = point.coord
        point_coord_PJCS = MsTrans.transPointCoord_CLCS2PJCS(point.coord, crs, M)  # from CLCS to PJCS
        Coordinates = [MetashapeType2numpyType([point_coord_CLCS[0], point_coord_CLCS[1], point_coord_CLCS[2]]),
                       MetashapeType2numpyType([point_coord_PJCS[0], point_coord_PJCS[1], point_coord_PJCS[2]])]
        Point[0] = Coordinates

        # Covariance
        if WriteCovariance:
            point_cov_CLCS = MetashapeType2numpyType(point.cov)  # return matrix with shape 3x3
            point_cov_LSECS = R * point.cov * R.t()  # # Transform point covariance matrix from CLCS to PJCS
            StandardDeviation = [1000 * math.sqrt(point_cov_LSECS[0, 0]),
                                 1000 * math.sqrt(point_cov_LSECS[1, 1]),
                                 1000 * math.sqrt(point_cov_LSECS[2, 2])]
            Point[1] = point_cov_CLCS
            Point[2] = MetashapeType2numpyType(StandardDeviation)  # return ndarray([sigX,sigY,sigZ])

        # Quality
        if WriteQuality:
            PointCriterion = PointsCriterion[point_id]
            Quality = [PointCriterion[0],
                       PointCriterion[1],
                       PointCriterion[2],
                       PointCriterion[3]]
            Point[3] = MetashapeType2numpyType(Quality)  # return ndarray([RE, RU, PA, IC])
        # State
        if point.selected:
            Point[4] = 1
        if point.valid:
            Point[5] = 1
        Point[6] = int(point.track_id)
        Points[point_id] = Point
    print('[Script][TimeCost]    DataImporting_Metashape.getComponent_Points:', time.perf_counter() - Functime, 'sec')
    return Points


def getComponent_Sensors(chunk):
    Functime = time.perf_counter()
    Sensors = {}
    for sensor_id, sensor in enumerate(chunk.sensors):
        calibration_params = MetashapeType2numpyType([sensor.calibration.f,
                                                      sensor.calibration.cx, sensor.calibration.cy,
                                                      sensor.calibration.b1, sensor.calibration.b2,
                                                      sensor.calibration.k1, sensor.calibration.k2,
                                                      sensor.calibration.k3, sensor.calibration.k4,
                                                      sensor.calibration.p1, sensor.calibration.p2,
                                                      sensor.calibration.p3, sensor.calibration.p4,
                                                      sensor.calibration.height, sensor.calibration.width])
        calibration_matrix = (sensor.calibration.covariance_matrix)
        Sensor = [calibration_params, calibration_matrix]
        Sensors[sensor_id] = Sensor
    print('[Script][TimeCost]    DataImporting_Metashape.getComponent_Sensors:', time.perf_counter() - Functime, 'sec')
    return Sensors


def getComponent_Projections(chunk):
    '''
    Projections = {identifier:[[u,v,size,track_id],...],...}
    projection_info = [projection.coord[0], projection.coord[1], projection.size].
    '''
    Functime = time.perf_counter()
    Projections = {}
    projections = chunk.point_cloud.projections
    for camera in chunk.cameras:
        if camera.type == Metashape.Camera.Type.Keyframe:
            print('[Script]    the sensor is keyframe!:', camera.label)
            continue  # skipping Keyframes
        elif not camera.transform:
            print('[Script]    camera has no transform!:', camera.label)
            continue  # skipping no transform
        elif not camera.enabled:
            print('[Script]    camera is not enabled!', camera.label)
            continue  # skipping camera not used
        # construct key
        identifier = MsCam.getCameraIdentifier(camera)
        # construct value
        Projection = []
        for projection in projections[camera]:
            projection_info = [projection.coord[0], projection.coord[1], projection.size, projection.track_id]
            Projection.append(projection_info)
        Projections[identifier] = Projection
    print('[Script][TimeCost]    DataImporting_Metashape.getComponent_Projections:', time.perf_counter() - Functime, 'sec')
    return Projections


def getComponent_Coordinates(chunk, PJCS_EPSG):
    '''
    Note: the default coordinate system of document should be set as PJCS.
    PJCS_EPSG is the authority of PJCS, such as 32650 (int or str)
    '''
    Functime = time.perf_counter()
    # Coordinate attribute
    crs_name = chunk.crs.name  # type string
    # Coordinate wkt format
    crs_wkt = chunk.crs.wkt  # PJCS should be the default
    if crs_wkt[:6] == 'PROJCS':
        crs_PJCS_wkt = crs_wkt
    else:
        crs_PJCS_wkt = Metashape.CoordinateSystem('EPSG::' + str(PJCS_EPSG)).wkt
    crs_GCCS_wkt = chunk.crs.geoccs.wkt
    crs_GDCS_wkt = chunk.crs.geogcs.wkt
    CoordinateAttribute = [crs_name, crs_PJCS_wkt, crs_GCCS_wkt, crs_GDCS_wkt]

    # Coordinate transform matrix
    # point in CLCS multiplied by this matrix could transform it to GCCS
    M = chunk.transform.matrix
    # point in GCCS multiplied by this matrix could transform it to LSECS
    LSE = chunk.crs.localframe(M.mulp(chunk.region.center))
    # point in CLCS multiplied by this matrix could transform it to LSECS
    T = LSE * M
    # rotation component of CLCS in LSECS
    if chunk.transform.scale:
        R = chunk.transform.scale * T.rotation()
    else:
        R = T.rotation()
    CoordinateTransform = [MetashapeType2numpyType(M),
                           MetashapeType2numpyType(LSE),
                           MetashapeType2numpyType(T),
                           MetashapeType2numpyType(R),
                           M.scale(),
                           MetashapeType2numpyType(M.inv())]
    print('[Script][TimeCost]    DataImporting_Metashape.getComponent_Coordinates:', time.perf_counter() - Functime, 'sec')
    return CoordinateTransform, CoordinateAttribute


def getCameraIdentifier(camera):
    '''
    make identifier
    '''
    # '2021:04:04 13:40:13' -> 20210404134013
    identifier = ''
    for s in camera.photo.meta['Exif/DateTimeOriginal']:
        if s.isdigit():
            identifier += s
    # 'DSC00009.JPG' -> 09  (# camera.label is DSC00004 in Metashape 1.8.5)
    # identifier += camera.label[-6:-4]
    identifier += camera.label[-2:]
    return int(identifier)


def getCameraPaths(chunk):
    '''
    Creat camera path dictionary for get identifier from camera path.
    this function could be used to link the images processed in other software to Metashape.chunk
    '''
    Functime = time.perf_counter()
    cameraPaths = {}
    for camera_id, camera in enumerate(chunk.cameras):
        # skipping Keyframes
        if camera.type == Metashape.Camera.Type.Keyframe:
            print('[Script]    the sensor is keyframe!:', camera.label)
            continue
        # make identifier
        identifier = getCameraIdentifier(camera)
        # write [key:value->identifier:camera_id] to camera_ids
        cameraPaths[camera.photo.path] = identifier
    print('[Script][TimeCost]    DataImporting_Metashape.getCameraPaths:', time.perf_counter() - Functime, 'sec')
    return cameraPaths


def getCameraIds(chunk):
    '''
    When creating several projects with subsets of the same image set, the same camera_id across projects may not
    indicate the same image. To handle this case, we creat a unique image identifier for each image using its
    photographing date and file size.

    camera_ids = {}:
        dict that return camera_id by given camera_Identifier
    '''
    Functime = time.perf_counter()
    camera_ids = {}
    for camera_id, camera in enumerate(chunk.cameras):
        # skipping Keyframes
        if camera.type == Metashape.Camera.Type.Keyframe:
            print('[Script]    the sensor is keyframe!:', camera.label)
            continue
        identifier = MsCam.getCameraIdentifier(camera)
        # write [key:value->identifier:camera_id] to camera_ids
        camera_ids[identifier] = camera_id
    print('[Script][TimeCost]    DataImporting_Metashape.getCameraIds:', time.perf_counter() - Functime, 'sec')
    return camera_ids


def getPointIds(chunk):
    '''
    Since Metashape provide attribute point.track_id, but couldn't access point_id by track_id,
    we creat a list that return point_id by given track_id.
    '''
    Functime = time.perf_counter()
    point_ids = [-1] * len(chunk.point_cloud.tracks)
    for point_id, point in enumerate(chunk.point_cloud.points):
        point_ids[point.track_id] = point_id
    print('[Script][TimeCost]    DataImporting_Metashape.getPointIds:', time.perf_counter() - Functime, 'sec')
    return point_ids


def getTracks(chunk):
    '''
    Tracks = [Track,Track,...,Track], a list that return Track by given track_id
        ,where Track = [view,view,...view],
        ,where view = [camera_id,projection_info]
        ,where projection_info = [projection.coord[0], projection.coord[1], projection.size, projection.track_id]
    '''
    Functime = time.perf_counter()
    projections = chunk.point_cloud.projections
    Tracks = [[] for i in range(len(chunk.point_cloud.tracks))]
    for camera_id, camera in enumerate(chunk.cameras):
        if camera.type == Metashape.Camera.Type.Keyframe:
            print('[Script]    the sensor is keyframe!:', camera.label)
            continue  # skipping Keyframes
        elif not camera.transform:
            print('[Script]    camera has no transform!:', camera.label)
            continue  # skipping no transform
        elif not camera.enabled:
            print('[Script]    camera is not enabled!', camera.label)
            continue  # skipping camera not used
        for projection in projections[camera]:
            view = [camera_id, [projection.coord[0], projection.coord[1], projection.size, projection.track_id]]
            Tracks[projection.track_id].append(view)
    print('[Script][TimeCost]    DataImporting_Metashape.getTracks:', time.perf_counter() - Functime, 'sec')
    return Tracks

def getTrackIds(chunk):
    '''
    creat a list that return track_id by given point_id.
    '''
    Functime = time.perf_counter()
    track_ids = [-1] * len(chunk.point_cloud.points)
    for point_id, point in enumerate(chunk.point_cloud.points):
        track_ids[point_id] = point.track_id
    print('[Script][TimeCost]    DataImporting_Metashape.getTrackIds:', time.perf_counter() - Functime, 'sec')
    return track_ids

def getValidTracks(chunk, point_ids, Tracks):
    '''
    validTracks is a part of Tracks, which contains only the Track that generate point successfully
    Considering the index of validTracks is meaningless, the validTracks_ids is used to record track_id of valid Track
    '''
    Functime = time.perf_counter()
    points = chunk.point_cloud.points
    validTracks = []
    validTracks_ids = []
    for track_id, point_id in enumerate(point_ids):
        if point_id == -1 or not points[point_id].valid:
            continue
        else:
            validTracks.append(Tracks[track_id])
            validTracks_ids.append(track_id)
    print('[Script][TimeCost]    DataImporting_Metashape.getValidTracks:', time.perf_counter() - Functime, 'sec')
    return validTracks, validTracks_ids


def getMatches(chunk, point_ids, Tracks):
    '''
    return:
        pair_ids = {}:
            dict that return [all, valid, invalid] by given (camera_id,camera_id);
        match_ids = {}:
            dict that return [valid matches, invalid matches] by given (camera_id,camera_id);
            where matches = [match,match,...,match],
            where match = (projection_info,projection_info)
            where projection_info = [projection.coord[0], projection.coord[1], projection.size].

    Note that the pair_ids and match_ids is completely the same with the matches viewed in GUI.
    To our knowledge, the matches here are only those used to construct Track, and the valid match means the Track
    triangulated successfully.
    So the valid or invalid flag only refers whether this match is USED, the process of geometric verified, misalign of
    Track and RANSAC based Track filtering might not considered here.
    '''
    Functime = time.perf_counter()
    cameras = chunk.cameras
    points = chunk.point_cloud.points

    # Creat full connection pairs
    pair_ids = {}
    match_ids = {}
    for camera_id1 in range(len(cameras) - 1):
        for camera_id2 in range(camera_id1 + 1, len(cameras)):
            pair_ids[(camera_id1, camera_id2)] = [0, 0, 0]  # [all, valid, invalid]
            match_ids[(camera_id1, camera_id2)] = [[], []]  # [valid matches, invalid matches]

    for track_id, Track in enumerate(Tracks):
        point_id = point_ids[track_id]
        for i, view1 in enumerate(Track[:-1]):
            for view2 in Track[i + 1:]:
                camera_id1 = view1[0]
                camera_id2 = view2[0]
                if camera_id1 < camera_id2:
                    # pair = (camera_epoch1,camera_epoch2)
                    pair = (camera_id1, camera_id2)
                    # match = ([proj_u,proj_v,size,track_id],[proj_u,proj_v,size,track_id])
                    match = (view1[1], view2[1])
                elif camera_id1 > camera_id2:
                    pair = (camera_id2, camera_id1)
                    match = (view2[1], view1[1])
                else:
                    raise Exception("[Script]    creating pair:a Track contains same camera!")

                # add matches number to dirt
                pair_ids[pair][0] += 1
                if point_id == -1 or not points[point_id].valid:
                    pair_ids[pair][2] += 1
                    match_ids[pair][1].append(match)
                else:
                    pair_ids[pair][1] += 1
                    match_ids[pair][0].append(match)
    print('[Script][TimeCost]    DataImporting_Metashape.getMatches:', time.perf_counter() - Functime, 'sec')
    return pair_ids, match_ids


