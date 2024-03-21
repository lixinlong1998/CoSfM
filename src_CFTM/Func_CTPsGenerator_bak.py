import math
import time
import multiprocessing
import numpy as np

import src_CFTM.Func_Geometry as Geom
import src_CFTM.Func_Camera as Cam
import src_CFTM.Func_Transform as Trans
import src_CFTM.Func_Match as Matc
import src_CFTM.Func_Track as Trac
import src_CFTM.Func_Files as Wrif
import src_CFTM.Func_CommonTiePoints as Comm

'''
step1: Creat grid to separate sparse cloud points.
step2: Construct cube for each cell of grid.
step3: Filter visible cameras by points in cell.
step4: Project cube to visible cameras to get specific sets of feature points.
step5: Match features cross all visible cameras, and validate the match if cameras belong to different epoch.
step6: Construct epoch tracks from matches.
step7: Filter the tracks by analysing number of CommonMatches and distance of features.
step8: Construct markers based on epoch Track.
v2中对cube生成的track包含了太多误匹配点,并且第一个track通常很大，而后续的track几乎很少有合适的。
v3还是选择v1的思路，但是添加滤波器。
'''


def GenerateCTPs(grid_id, corners, Camera_IdList, Cameras, camera_ids, CamerasEpoch, Sensors,
                 CoordinateTransform, CoordinateAttribute, chunkKeypoints, chunkDescriptors, path, arguments):
    # get arguments
    criterias1 = arguments[0]
    criterias2 = arguments[1]
    threshold_RepError2 = arguments[2]
    ratio_inlier_init = arguments[3]
    confidence = arguments[4]
    threshold_Distance = arguments[5] / CoordinateTransform[4]

    # start
    starttime0 = time.perf_counter()
    print('[Script]    [{}]Generating CTPs:'.format(grid_id), grid_id)

    # Step1: project cube to visible cameras to extract feature subset
    print('[Script]    [{}]    Extract cameras features...'.format(grid_id))
    starttime1 = time.perf_counter()
    features_cube_e1 = []
    features_cube_e2 = []
    for camera_id in Camera_IdList:
        Camera = Cameras[camera_id]
        identifier = Camera[-1]
        # project cube to each visible camera to get polygon of cube
        cube2DVertexes = projectCube(Camera, corners, Sensors, CoordinateTransform, CoordinateAttribute)
        polygon = Geom.encloseCubePolygon(cube2DVertexes)
        # extract keypoints by polygon of cube
        features_cube_camera_id = extractFeatures_byMask(identifier, polygon, chunkKeypoints)
        # [Exception Handling]: In case there is no feature located in the projection polygon
        if len(features_cube_camera_id) == 1:
            # skip this camera
            continue
        # append subset of features to corresponding epoch list
        if CamerasEpoch[camera_id] == 0:
            features_cube_e1.append(features_cube_camera_id)
        elif CamerasEpoch[camera_id] == 1:
            features_cube_e2.append(features_cube_camera_id)
        else:
            raise Exception("[Script]    this script could only deal with pairwise co-align!")
    print('[Script]    [{}]    cameras number in e1:'.format(grid_id), len(features_cube_e1))
    print('[Script]    [{}]    cameras number in e2:'.format(grid_id), len(features_cube_e2))
    # [Exception Handling]: In case there only one epoch cameras could observe this cube
    if len(features_cube_e1) < 2 or len(features_cube_e2) < 2:
        '''
        'or' or 'and' indicate two different robust degrees;
          'or' means no matter how many cameras are visible in one epoch, as long as the other epoch has
                     only one camera visible, this cube will be ignored;
          'and' means when each epoch only contained one visible camera, this cube then will be ignored.
        '''
        print('[Script]    [{}]Too less cameras could observe this cube!'.format(grid_id))
        print('[Script]    [{}]Skip grid!'.format(grid_id))
        return None
    assert len(features_cube_e1) >= 2
    assert len(features_cube_e2) >= 2
    print('[Script]    [{}]    [Script][TimeCost]    :'.format(grid_id), time.perf_counter() - starttime1)

    # Step2: Match features in each epoch
    starttime1 = time.perf_counter()
    print('[Script]    [{}]    Matching features...'.format(grid_id))
    matches_cube_e1 = Matc.innerMatch(features_cube_e1, chunkKeypoints, chunkDescriptors, Cameras, Sensors, camera_ids,
                                      criterias1)
    matches_cube_e2 = Matc.innerMatch(features_cube_e2, chunkKeypoints, chunkDescriptors, Cameras, Sensors, camera_ids,
                                      criterias1)
    print('[Script]    [{}]    Number of camera pairs in e1:'.format(grid_id), len(matches_cube_e1))
    print('[Script]    [{}]    Number of camera pairs in e2:'.format(grid_id), len(matches_cube_e2))
    # [Exception Handling]: In case there is no valid match left in matches_cube_epoch
    if len(matches_cube_e1) == 0 or len(matches_cube_e2) == 0:
        print('[Script]    [{}]This cube has no matched camera pair!'.format(grid_id))
        print('[Script]    [{}]Skip grid:', grid_id)
        return None
    # [export]

    print('[Script]    [{}]    [Script][TimeCost]    :'.format(grid_id), time.perf_counter() - starttime1)

    # Step3: Construct tracks from matches for each epoch
    starttime1 = time.perf_counter()
    print('[Script]    [{}]    Constructing tracks...'.format(grid_id))
    Tracks_Generated_e1, Tracks_Matches_e1 = Trac.generateTracks(matches_cube_e1)
    Tracks_Generated_e2, Tracks_Matches_e2 = Trac.generateTracks(matches_cube_e2)
    Tracks_length1 = [len(track_i) for track_i in Tracks_Generated_e1]
    Tracks_length2 = [len(track_i) for track_i in Tracks_Generated_e2]
    print('[Script]    [{}]    tracks in e1:'.format(grid_id), max(Tracks_length1), np.mean(Tracks_length1),
          min(Tracks_length1))
    print('[Script]    [{}]    tracks in e2:'.format(grid_id), max(Tracks_length2), np.mean(Tracks_length2),
          min(Tracks_length2))
    print('[Script]    [{}]    Number of generated tracks in e1:'.format(grid_id), len(Tracks_Generated_e1))
    print('[Script]    [{}]    Number of generated tracks in e2:'.format(grid_id), len(Tracks_Generated_e2))
    # [Exception Handling]: In case there is no tracks generated
    if len(Tracks_Generated_e1) == 0 or len(Tracks_Generated_e2) == 0:
        print('[Script]    [{}]This cube has no tracks generated!'.format(grid_id))
        print('[Script]    [{}]Skip grid!'.format(grid_id))
        return None
    print('[Script]    [{}]    [Script][TimeCost]    :'.format(grid_id), time.perf_counter() - starttime1)

    # Step4: Filter Tracks
    starttime1 = time.perf_counter()
    print('[Script]    [{}]    Filtering tracks...'.format(grid_id))
    Tracks_Filtered_e1, Tracks_Quality_e1 = Trac.filterTracks(grid_id, Tracks_Generated_e1, Cameras, Sensors,
                                                              chunkKeypoints, ratio_inlier_init, confidence,
                                                              criterias2, threshold_RepError2)
    Tracks_Filtered_e2, Tracks_Quality_e2 = Trac.filterTracks(grid_id, Tracks_Generated_e2, Cameras, Sensors,
                                                              chunkKeypoints, ratio_inlier_init, confidence,
                                                              criterias2, threshold_RepError2)
    print('[Script]    [{}]    Number of tracks left in e1:'.format(grid_id), len(Tracks_Filtered_e1))
    print('[Script]    [{}]    Number of tracks left in e2:'.format(grid_id), len(Tracks_Filtered_e2))
    # [Exception Handling]: In case there is no epoch Track left
    if len(Tracks_Filtered_e1) == 0 or len(Tracks_Filtered_e2) == 0:
        print('[Script]    [{}]This cube has no epoch Track left after filter!'.format(grid_id))
        print('[Script]    [{}]Skip grid!'.format(grid_id))
        return None
    print('[Script]    [{}]    [Script][TimeCost]    :'.format(grid_id), time.perf_counter() - starttime1)

    # Step5: Match tracks from different epochs with robust filters
    starttime1 = time.perf_counter()
    print('[Script]    [{}]    Matching tracks...'.format(grid_id))
    Tracks_Matched, Tracks_Report = Trac.matchTracks(Tracks_Filtered_e1, Tracks_Filtered_e2,
                                                     Tracks_Quality_e1, Tracks_Quality_e2,
                                                     Cameras, Sensors, chunkDescriptors, threshold_Distance,crossCheck=True)
    print('[Script]    [{}]    Number of epoch tracks:'.format(grid_id), len(Tracks_Matched))
    # [Exception Handling]: In case there is no epoch Track left
    if len(Tracks_Matched) == 0:
        print('[Script]    [{}]This cube has no epoch Track!'.format(grid_id))
        print('[Script]    [{}]Skip grid!'.format(grid_id))
        return None
    print('[Script]    [{}]    [Script][TimeCost]    :'.format(grid_id), time.perf_counter() - starttime1)

    # Step6: Write Tracks_Best to file.txt with form of Metashape.Markers
    starttime1 = time.perf_counter()
    print('[Script]    [{}]    Writing tracks...'.format(grid_id))
    Wrif.writeCommonTracks(path, Tracks_Matched, Tracks_Report)
    print('[Script]    [{}]    Number of Marker Points of this cube:'.format(grid_id), len(Tracks_Matched))
    print('[Script]    [{}]    [Script][TimeCost]    :'.format(grid_id), time.perf_counter() - starttime1)

    print('[Script]    [{}][Script][TimeCost]    :'.format(grid_id), time.perf_counter() - starttime0)
    return None


def buildMarkersGrid(Boundary_PJCS, girdsize, MarkerPoints, MarkerTracks, MarkerInforms):
    '''
    creat gird
    input:
        MarkerPoints = [MarkerPoint,MarkerPoint,...,MarkerPoint]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
            ,where Point_Coord_CLCS = list[X,Y,Z] in CLCS
        MarkerTracks = [track,track,...,track]
            ,where track = [view,view,...,view]
            ,where view = [camera_id, projection_info, *]
            ,where projection_info = [u, v, *]
            ,where * = [identifier, keypoints_id]
        MarkerInforms = [Information,Information,...,Information]
            ,where Information=[[point_id,Track_id],] or [[similarity, angles, RepErrors],]
    output:
        MarkersGrid = {(r, c):[Marker,...]}
            ,where Marker = [MarkerPoint,track,Information]
    '''
    # get boundary
    boundary_X_MAX, boundary_X_MIN, boundary_Y_MAX, boundary_Y_MIN = Boundary_PJCS
    # creat empty PointsGrid and ICTPsGrid
    rows = math.ceil((boundary_Y_MAX - boundary_Y_MIN) / girdsize)
    cols = math.ceil((boundary_X_MAX - boundary_X_MIN) / girdsize)
    print('[Script]    rows, cols:', rows, cols)
    # generate points dictionary corresponds to girds
    MarkersGrid = {}
    for r in range(rows):
        for c in range(cols):
            MarkersGrid[(r, c)] = []
    # append point_id to correspond each cell of girds
    for i, MarkerPoint in enumerate(MarkerPoints):
        Point_Coord_PJCS = MarkerPoint[1]
        Coordinates_PJCS_X = Point_Coord_PJCS[0]
        Coordinates_PJCS_Y = Point_Coord_PJCS[1]
        r = int((boundary_Y_MAX - float(Coordinates_PJCS_Y)) / girdsize)  # row distance to origin(top-left)
        c = int((float(Coordinates_PJCS_X) - boundary_X_MIN) / girdsize)  # col distance to origin(top-left)
        # combine a marker object
        Marker = [MarkerPoint, MarkerTracks[i], MarkerInforms[i]]
        MarkersGrid[(r, c)].append(Marker)
    return MarkersGrid


def buildPointsGrid(Boundary_PJCS, Points, ICTPsCoord_PJCS, girdsize):
    '''
    creat gird
    '''
    # get boundary
    boundary_X_MAX, boundary_X_MIN, boundary_Y_MAX, boundary_Y_MIN = Boundary_PJCS
    # creat empty PointsGrid and ICTPsGrid
    rows = math.ceil((boundary_Y_MAX - boundary_Y_MIN) / girdsize)
    cols = math.ceil((boundary_X_MAX - boundary_X_MIN) / girdsize)
    print('[Script]    rows, cols:', rows, cols)
    ICTPsGrid = np.mat(np.zeros([rows, cols], dtype=np.uintc, order='C'))
    # generate points dictionary corresponds to girds
    PointsGrid = {}
    for r in range(rows):
        for c in range(cols):
            PointsGrid[(r, c)] = []
    # append point_id to correspond each cell of girds
    for point_id, Point in Points.items():
        Coordinates_PJCS = Point[0][1]
        Coordinates_PJCS_X = Coordinates_PJCS[0]
        Coordinates_PJCS_Y = Coordinates_PJCS[1]
        r = int((boundary_Y_MAX - float(Coordinates_PJCS_Y)) / girdsize)  # row distance to origin(top-left)
        c = int((float(Coordinates_PJCS_X) - boundary_X_MIN) / girdsize)  # col distance to origin(top-left)
        PointsGrid[(r, c)].append(point_id)
    # judge whether cell contains CTP
    for MCTPCoordXY in ICTPsCoord_PJCS:
        MCTPCoord_X = MCTPCoordXY[0]
        MCTPCoord_Y = MCTPCoordXY[1]
        r = int((boundary_Y_MAX - float(MCTPCoord_Y)) / girdsize)  # row distance to origin(top-left)
        c = int((float(MCTPCoord_X) - boundary_X_MIN) / girdsize)  # col distance to origin(top-left)
        ICTPsGrid[r, c] += 1
    return PointsGrid, ICTPsGrid


def buildPointsCubes(Boundary_PJCS, Points, TracksEpoch, PointsGrid_NoICTPs, girdsize):
    '''
    Creat cubes for those cells contains points from different epochs.
    The upper and lower plane of cube is determined by the maximum and minimum of Z value with double std.
    '''
    # extract point Z coordinate and its standard deviation from Points
    PointsZ = {}
    PointsEpoch = {}
    for point_id, Point in Points.items():
        PointsZ[point_id] = [Point[0][1][2], Point[2][2]]  # Z(m),sigZ(mm)
        track_id = Point[6]
        signal = Comm.analyseSignal_Epoch(TracksEpoch[track_id])
        PointsEpoch[point_id] = signal

    # Creat cubes
    Cubes = {}
    for grid_id, Point_IdList in PointsGrid_NoICTPs.items():
        # skip cell without points
        if len(Point_IdList) == 0:
            continue
        # skip cell if only contains points from one epoch(points number below 10 will be see as no points)
        signals = [PointsEpoch[point_id] for point_id in Point_IdList]
        epochNum = np.sum(np.asarray(signals), axis=0)
        signal=[]
        for i in epochNum:
            if i <= 10:
                signal.append(0)
            else:
                signal.append(1)
        if sum(signal) == 1:
            continue
        # get coordinates of four conners of cell
        r, c = grid_id
        boundary_X_MAX, boundary_X_MIN, boundary_Y_MAX, boundary_Y_MIN = Boundary_PJCS
        topLeft = [boundary_X_MIN + c * girdsize, boundary_Y_MAX - r * girdsize]
        topRight = [boundary_X_MIN + (c + 1) * girdsize, boundary_Y_MAX - r * girdsize]
        leftBottom = [boundary_X_MIN + c * girdsize, boundary_Y_MAX - (r + 1) * girdsize]
        rightBottom = [boundary_X_MIN + (c + 1) * girdsize, boundary_Y_MAX - (r + 1) * girdsize]
        # get the upper and lower coordinate of cube
        Cube_Z = []
        for point_id in Point_IdList:
            Z, sigmaZ = PointsZ[point_id]  # in (m) and (mm) respectively
            Cube_Z.append(Z - 3 * sigmaZ / 1000)
            Cube_Z.append(Z + 3 * sigmaZ / 1000)
        upper = max(Cube_Z)
        lower = min(Cube_Z)

        Cubes[grid_id] = [topLeft, topRight, leftBottom, rightBottom, upper, lower]
    print('cubes number:', len(Cubes))
    return Cubes


def getGrids_NoICTPs(PointsGrid, ICTPsGrid):
    PointsGrid_NoICTPs = {}
    rows, cols = ICTPsGrid.shape
    for r in range(rows):
        for c in range(cols):
            if ICTPsGrid[r, c] == 0:
                PointsGrid_NoICTPs[(r, c)] = PointsGrid[(r, c)]
    return PointsGrid_NoICTPs


def getGridDots(boundarys, girdsize):
    # get boundary
    boundary_X_MAX, boundary_X_MIN, boundary_Y_MAX, boundary_Y_MIN = boundarys
    # creat empty pointsGrid and ICTPsGrid
    rows = math.ceil((boundary_Y_MAX - boundary_Y_MIN) / girdsize)
    cols = math.ceil((boundary_X_MAX - boundary_X_MIN) / girdsize)
    GridDots = []
    for r in range(rows + 1):
        for c in range(cols + 1):
            GridDots.append([boundary_X_MIN + c * girdsize, boundary_Y_MAX - r * girdsize])
    return GridDots


def projectCube(Camera, corners, Sensors, CoordinateTransform, CoordinateAttribute):
    '''
    input:
        corners = [topLeft, topRight, leftBottom, rightBottom, upper, lower]
        Camera=[transform, reference, covariance, sensors_id, state,path, identifier]
    output:
        vertexes = np.2darray([[X,Y],[X,Y],[X,Y],[X,Y],[X,Y],[X,Y],[X,Y],[X,Y]])
        ,where the points are unordered.
    '''
    UpTopLeft = [corners[0][0], corners[0][1], corners[4]]
    UpTopRight = [corners[1][0], corners[1][1], corners[4]]
    UpLeftBottom = [corners[2][0], corners[2][1], corners[4]]
    UpRightBottom = [corners[3][0], corners[3][1], corners[4]]
    DownTopLeft = [corners[0][0], corners[0][1], corners[5]]
    DownTopRight = [corners[1][0], corners[1][1], corners[5]]
    DownLeftBottom = [corners[2][0], corners[2][1], corners[5]]
    DownRightBottom = [corners[3][0], corners[3][1], corners[5]]
    cube3DVertexes = [UpTopLeft, UpTopRight, UpLeftBottom, UpRightBottom,
                      DownTopLeft, DownTopRight, DownLeftBottom, DownRightBottom]
    cube2DVertexes = []
    for Vertex3D_PJCS in cube3DVertexes:
        Vertex3D_CLCS = Trans.transPointCoord_PJCS2CLCS(Vertex3D_PJCS, CoordinateTransform, CoordinateAttribute)
        Vertex2D = Cam.cameraProjector(Vertex3D_CLCS, Camera, Sensors)
        cube2DVertexes.append([Vertex2D[0], Vertex2D[1]])
    return np.asarray(cube2DVertexes)


def extractFeatures_byMask(identifier, polygon, chunkKeypoints):
    '''
    extract keypoints by polygon of cube
    input:
        chunkKeypoints = {identifier:2darray.shape(numKP,6).dtype(float32),...}
        polygon = np.2darray([[X,Y],[X,Y],...,[X,Y]])
        identifier is unique number for camera in document
    output:
        features_cube_camera_id = [identifier,keypoint_id,keypoint_id,...,keypoint_id]
    '''
    features_cube_camera_id = [identifier]
    for keypoint_id, keypoint in enumerate(chunkKeypoints[identifier]):
        # chunkKeypoints[identifier] : 2darray.shape(numKP,6).dtype(float32)
        keypoint_coord = keypoint[:2]  # 1darray.shape(2,)
        if not Geom.isInPolygon(keypoint_coord, polygon):
            continue
        features_cube_camera_id.append(keypoint_id)
    return features_cube_camera_id
