from config import *
import Metashape
import csv
import numpy as np
from src_CFTM import ConnectData_Metashape
from src_CFTM import Func_Statistic as FcSta
from src_CFTM import Func_Files as FcFile
import src_Metashape.FuncMs_Triangulate as MsTri
import src_Metashape.FuncMs_Transform as MsTrans

'''
Metashape.Markers
如果想要直接分析Markers，使用getMarkers，exportMarkersData_Analyse，importMarkersData_Analyse
如果想要添加额外的信息到Markers中，使用getMarkers，exportMarkersData_Add，importMarkersData_Add
'''

'''添加Markers'''


# 从chunk中的Metashape.Markers读取其【数据添加格式】
def getMarkersData_Add(chunk, MarkerGroupName, mode='both'):
    '''
    get track, point and information from Metashape.Markers in document
    input:
        Metashape.Document.chunk.markers
        MarkerGroupName
        mode = 'check', 'uncheck', 'both'
    output:
        MarkerPoints = [MarkerPoint,MarkerPoint,...,MarkerPoint]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
            ,where Point_Coord_CLCS = list[X,Y,Z] in CLCS
        MarkerTracks = {track,track,...,track}
            ,where Track = [view,view,...,view]
            ,where view = [camera_id, projection_info]
            ,where projection_info = [u, v]
        MarkerInforms = [MarkerInform,MarkerInform,...,MarkerInform]
            ,where MarkerInform = [marker.label, marker.key]
    '''
    if MarkerGroupName:
        for markerGroup in chunk.marker_groups:
            if markerGroup.label != MarkerGroupName:
                continue
            # Now, we get the markerGroup of given MarkerGroupName
            MarkerPoints = {}
            MarkerTracks = {}
            MarkerInforms = {}
            for marker_id, marker in enumerate(chunk.markers):
                if marker.group != markerGroup:
                    continue
                if not MarkerMode(marker, mode):
                    continue
                Point_Coord_CLCS = marker.position
                Point_Coord_PJCS = marker.reference.location
                MarkerPoints[marker_id] = [[Point_Coord_CLCS[0], Point_Coord_CLCS[1], Point_Coord_CLCS[2]],
                                           [Point_Coord_PJCS[0], Point_Coord_PJCS[1], Point_Coord_PJCS[2]]]
                MarkerTracks[marker_id] = getMarkerTrack(marker, chunk.cameras, Index='camera_id')
                MarkerInforms[marker_id] = [marker.label, marker.key]
    else:
        MarkerPoints = {}
        MarkerTracks = {}
        MarkerInforms = {}
        for marker_id, marker in enumerate(chunk.markers):
            if not MarkerMode(marker, mode):
                continue
            Point_Coord_CLCS = marker.position
            Point_Coord_PJCS = marker.reference.location
            MarkerPoints[marker_id] = [[Point_Coord_CLCS[0], Point_Coord_CLCS[1], Point_Coord_CLCS[2]],
                                       [Point_Coord_PJCS[0], Point_Coord_PJCS[1], Point_Coord_PJCS[2]]]
            MarkerTracks[marker_id] = getMarkerTrack(marker, chunk.cameras, Index='camera_id')
            MarkerInforms[marker_id] = [marker.label, marker.key]
    return MarkerPoints, MarkerTracks, MarkerInforms


# 导出Metashape.Markers的【数据添加格式】
def exportMarkersData_Add(path, MarkerPoints, MarkerTracks, MarkerInforms):
    '''
    write track, point and information of Metashape.Marker to file.
    input:
        MarkerPoints = [MarkerPoint,MarkerPoint,...,MarkerPoint]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
        MarkerTracks = [MarkerTrack,MarkerTrack,...,MarkerTrack]
            ,where MarkerTrack = [view,view,...,view]
            ,where view = [camera_id, projection_info]
            ,where projection_info = [u, v]
        MarkerInforms = [MarkerInform,MarkerInform,...,MarkerInform]
            ,where MarkerInform = []
    output:
        .../markers.txt
    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    for i in range(len(MarkerTracks)):
        MarkerPoint = MarkerPoints[i]
        MarkerTrack = MarkerTracks[i]
        MarkerInform = MarkerInforms[i]
        fwriter.writerow([MarkerPoint, MarkerTrack, MarkerInform])
    File.close()


# 导入Metashape.Markers的【数据添加格式】
def importMarkersData_Add(path):
    '''
    read track, point and information of Metashape.Marker from file.
    input:
        .../markers.txt

    output:
        MarkerPoints = [MarkerPoint,MarkerPoint,...,MarkerPoint]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
        MarkerTracks = [MarkerTrack,MarkerTrack,...,MarkerTrack]
            ,where MarkerTrack = [view,view,...,view]
            ,where view = [camera_id, projection_info]
            ,where projection_info = [u, v]
        MarkerInforms = [MarkerInform,MarkerInform,...,MarkerInform]
            ,where MarkerInform = []
    '''
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        MarkerPoints = []
        MarkerTracks = []
        MarkerInforms = []
        for row in freader:
            # point coordinate
            MarkerPoints.append(eval(row[0]))
            # track
            MarkerTracks.append(eval(row[1]))
            # Informations
            MarkerInforms.append(eval(row[2]))
    return MarkerPoints, MarkerTracks, MarkerInforms


# 基于【数据添加格式】将Metashape.Markers添加到chunk中
def addMarkers(chunk, MarkerPoints, MarkerTracks, MarkerGroupName):
    '''
    add points into chunk with Markers format.
    input:
        chunk is Metashape.Document.chunk
        MarkerPoints = [MarkerPoint,MarkerPoint,...,MarkerPoint]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
            ,where Point_Coord_CLCS = list[X,Y,Z] in CLCS
        MarkerTracks = [track,track,...,track]
            ,where track = [view,view,...,view]
            ,where view = [camera_id, projection_info, *]
            ,where projection_info = [u, v, *]
            ,where * = [identifier, keypoints_id]
    output:
        Metashape.Document.chunk.markers
    '''
    cameras = chunk.cameras
    # creat new marker group
    MarkerGroup = chunk.addMarkerGroup()
    MarkerGroup.label = MarkerGroupName
    # creat new markers
    for markerName, MarkerPoint in enumerate(MarkerPoints):
        Point_Coord_CLCS = MarkerPoint[0]
        point_coord_CLCS = Metashape.Vector([Point_Coord_CLCS[0], Point_Coord_CLCS[1], Point_Coord_CLCS[2]])
        # creat new marker
        marker = chunk.addMarker(point_coord_CLCS, visibility=False)
        # assign marker's attributions
        marker.group = MarkerGroup
        marker.label = 'Marker_{0}'.format(markerName)
        # marker.reference.location = point_coord_PJCS
        marker.reference.accuracy = Metashape.Vector([10, 10, 10])  # in metres
        marker.reference.enabled = False
        # modify or append projection information, at the meanwhile, record the camera index list of MarkerTracks
        MarkerTracksCamera_IdList = []
        for view in MarkerTracks[markerName]:
            camera = cameras[view[0]]
            projection_info = view[1]
            projection_coord = Metashape.Vector([projection_info[0], projection_info[1]])
            # replace or append projection information of marker by MarkerTracks information
            marker.projections[camera] = Metashape.Marker.Projection(projection_coord, True)
            # record the camera index list of MarkerTracks
            MarkerTracksCamera_IdList.append(view[0])
        # delete projection of marker, which is not in MarkerTracks
        for camera, projection in marker.projections.items():
            if cameras.index(camera) not in MarkerTracksCamera_IdList:
                del marker.projections[camera]


# 将chunk中的点对象转换为Metashape.Markers的【数据添加格式】
def convertPoint2Marker(Points, Tracks, Points_IdList):
    '''
    input:
        chunk is Metashape.Document.chunk
    output:
        MarkerPoints = [MarkerPoint,MarkerPoint,...,MarkerPoint]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
            ,where Point_Coord_CLCS = list[X,Y,Z] in CLCS
        MarkerTracks = [track,track,...,track]
            ,where track = [view,view,...,view]
            ,where view = [camera_id, projection_info, *]
            ,where projection_info = [u, v, *]
            ,where * = [identifier, keypoints_id]
        MarkerInforms = [Information,Information,...,Information]
            ,where Information=[[point_id,Track_id],]

    '''
    MarkerPoints = []
    MarkerTracks = []
    MarkerInforms = []
    for point_id in Points_IdList:
        Point = Points[point_id]
        Track_id = Point[-1]
        Track = Tracks[Track_id]
        Point_Coord_CLCS = Point[0][0]
        Point_Coord_PJCS = Point[0][1]
        MarkerPoints.append([[Point_Coord_CLCS[0], Point_Coord_CLCS[1], Point_Coord_CLCS[2]],
                             [Point_Coord_PJCS[0], Point_Coord_PJCS[1], Point_Coord_PJCS[2]]])
        MarkerTracks.append(Track)
        MarkerInforms.append([[point_id, Track_id]])
    return MarkerPoints, MarkerTracks, MarkerInforms


# 将chunk中的track对象转换为Metashape.Markers的【数据添加格式】
def convertCommonTrack2Marker(CommonTracks, CommonTracksMatches, chunk):
    '''
    input:
        CommonTracks = [CommonTrack,CommonTrack,...,CommonTrack]
            ,where CommonTrack = [view,view,...,view]
            ,where view = [camera_id, projection_info, index_info]
            ,where projection_info = [u, v]
            ,where index_info = [identifier, keypoints_id, originTrack_id]
        CommonTracksMatches = [CommonTrackMatches,CommonTrackMatches,...,CommonTrackMatches]
            ,where CommonTrackMatches = [CrossMatch,…]
            ,where CrossMatch = [feature1, feature2, MatchQuality]
            ,where feature = [camera_id, keypoint_id]
            ,where MatchQuality = [distance,angles],distance is the L2 normal distance between descriptors of match.
    output:
        MarkerPoints = [MarkerPoint,MarkerPoint,...,MarkerPoint]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
            ,where Point_Coord_CLCS = list[X,Y,Z] in CLCS
        MarkerTracks = [MarkerTrack,MarkerTrack,...,MarkerTrack]
            ,where MarkerTrack = [view,view,...,view]
            ,where view = [camera_id, projection_info]
            ,where projection_info = [u, v]
        MarkerInforms = [Information,Information,...,Information]
            ,where Information=[[CrossMatch,…],]
    '''
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)
    cameras = chunk.cameras

    MarkerPoints = []
    MarkerTracks = []
    MarkerInforms = []
    for i, CommonTrack in enumerate(CommonTracks):
        Point_Coord_CLCS = MsTri.Triangulate(CommonTrack, cameras)
        Point_Coord_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point_Coord_CLCS, crs, M)
        MarkerPoints.append([[Point_Coord_CLCS[0], Point_Coord_CLCS[1], Point_Coord_CLCS[2]],
                             [Point_Coord_PJCS[0], Point_Coord_PJCS[1], Point_Coord_PJCS[2]]])
        MarkerTracks.append(CommonTrack)  # list[list]
        MarkerInforms.append([CommonTracksMatches[i]])
    return MarkerPoints, MarkerTracks, MarkerInforms


'''分析Markers'''


# 从chunk中的Metashape.Markers读取其【数据分析格式】
def getMarkersData_Analyse(chunk, MarkerGroupName, mode='both'):
    '''
    get track, point and information from Metashape.Markers in document
    input:
        Metashape.Document.chunk.markers
        camera_ids = {identifier:camera_id}
        MarkerGroupName
        mode = 'check', 'uncheck', 'both'
    output:
        MarkersDataAnalyse = {MarkerData,MarkerData,...,MarkerData}
            ,where MarkerData = [marker_id, marker.label, MarkerTrack]
            ,where MarkerTrack = [view,view,...,view]
            ,where view = [identifier, projection_info]
            ,where projection_info = [u, v]
    '''
    if MarkerGroupName:
        for markerGroup in chunk.marker_groups:
            if markerGroup.label != MarkerGroupName:
                continue
            # Now, we get the markerGroup of given MarkerGroupName
            MarkersDataAnalyse = []
            for marker_id, marker in enumerate(chunk.markers):
                if marker.group != markerGroup:
                    continue
                if not MarkerMode(marker, mode):
                    continue
                MarkersDataAnalyse.append(
                    [marker_id, marker.label, getMarkerTrack(marker, chunk.cameras, Index='identifier')])
    else:
        MarkersDataAnalyse = []
        for marker_id, marker in enumerate(chunk.markers):
            if not MarkerMode(marker, mode):
                continue
            MarkersDataAnalyse.append(
                [marker_id, marker.label, getMarkerTrack(marker, chunk.cameras, Index='identifier')])
    return MarkersDataAnalyse


# 从chunk中的Metashape.Markers导出【数据分析格式】
def exportMarkersData_Analyse(chunk, path, MarkerGroupName, mode='both'):
    '''
    export markers directly to file for recording basic data of marker.
    output:
        Markers = [Marker,Marker,...,Marker]
            ,where Marker = [marker_id, marker_label, MarkerTrack]
            ,where MarkerTrack = [view,view,...,view]
            ,where view = [identifier, projection_info]
            ,where projection_info = [u, v]
    '''
    starttime = time.perf_counter()
    # creat a file for writing results
    reportFile = open(path, "w")
    fwriter = csv.writer(reportFile, delimiter='\t', lineterminator='\n')
    if MarkerGroupName:
        for markerGroup in chunk.marker_groups:
            if markerGroup.label != MarkerGroupName:
                continue
            # Now, we get the markerGroup of given MarkerGroupName
            for marker_id, marker in enumerate(chunk.markers):
                if marker.group != markerGroup:
                    continue
                if not MarkerMode(marker, mode):
                    continue
                fwriter.writerow([marker_id, marker.label, getMarkerTrack(marker, chunk.cameras, Index='identifier')])
    else:
        for marker_id, marker in enumerate(chunk.markers):
            if not MarkerMode(marker, mode):
                continue
            fwriter.writerow([marker_id, marker.label, getMarkerTrack(marker, chunk.cameras, Index='identifier')])
    reportFile.close()
    print('[Script][TimeCost]    export CheckPoints:', time.perf_counter() - starttime)


# 导入Metashape.Markers的【数据分析格式】，以用于质量分析
def importMarkersData_Analyse(path, MarkerList=[]):
    '''
    import basic data of markers from file.
    output:
        Markers = [Marker,Marker,...,Marker]
            ,where Marker = [marker_id, marker_label, MarkerTrack]
            ,where MarkerTrack = [view,view,...,view]
            ,where view = [identifier, projection_info]
            ,where projection_info = [u, v]
    '''
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        Markers = []
        for row in freader:
            Marker_id = int(row[0])
            Marker_label = row[1]
            Marker_track = eval(row[2])
            # [marker_id, marker.label, getMarkerTrack_Analyse(marker)]
            Marker = [Marker_id, Marker_label, Marker_track]

            if len(MarkerList) == 0:
                # print(Marker_id, Marker_label)
                Markers.append(Marker)
            else:
                if Marker_id in MarkerList:
                    # print(Marker_id, Marker_label)
                    Markers.append(Marker)
    return Markers


# 基于Metashape.Markers的【数据分析格式】分析markers的质量
def getMarkersQua_Analyse(chunk, Markers, camera_ids, CamerasEpoch):
    '''
    input:
        Markers = [Marker,Marker,...,Marker]
            ,where Marker = [marker_id, marker_label, MarkerTrack]
            ,where MarkerTrack = [view,view,...,view]
            ,where view = [identifier, projection_info]
            ,where projection_info = [u, v]
        camera_ids = {identifier:camera_id}
    output:
        MarkersQua = [MarkerQua,MarkerQua,...,MarkerQua]
            ,where MarkerQua = [marker_id,
                                    0
                                Point3D_e1[0], Point3D_e1[1], Point3D_e1[2], Point3D_e2[0], Point3D_e2[1], Point3D_e2[2],
                                      1              2              3             4              5               6
                                residual_X, residual_Y, residual_Z, residual_P, residual_T,
                                     7          8           9           10          11
                                RepErr_v2_e1, RepErr_v3_e1, RepErr_v2_e2, RepErr_v3_e2,
                                     12            13            14            15
                                TriAng_Max, TriAng_Avg, TriAng_Min, ViewNum1, ViewNum2, EpoMatch_Num]
                                     16         17          18         19        20          21
    '''
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)
    cameras = chunk.cameras

    # calculate coordinate
    # here, we use CPs to represent an instance of Marker
    MarkersQua = []
    for marker_id, marker_label, MarkerTrack in Markers:
        # print('[Script]    Marker:', marker_id)

        # split MarkerTrack by epoch and convert identifier to camera_id
        Track1 = []
        Track2 = []
        for view in MarkerTrack:
            # try to read the camera of this view, if the camera is not exist, this view will be skipped.
            camera_id = camera_ids[int(view[0])]
            projection_info = view[1]
            if CamerasEpoch[camera_id] == 0:
                Track1.append([camera_id, projection_info])
            elif CamerasEpoch[camera_id] == 1:
                Track2.append([camera_id, projection_info])
            else:
                raise Exception("[Script]    the given epochs are not include all dates")

        # skip the case that sub-Track not contains enough images
        ViewNum1 = len(Track1)
        ViewNum2 = len(Track2)
        if ViewNum1 < 2 or ViewNum2 < 2:
            print('[Script]    Marker:', marker_id, '  skipped')
            continue

        # Triangulate
        Point3D_e1_CLCS = MsTri.Triangulate(Track1, cameras)
        Point3D_e2_CLCS = MsTri.Triangulate(Track2, cameras)
        Point3D_e1_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e1_CLCS, crs, M)
        Point3D_e2_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e2_CLCS, crs, M)

        # analyse quality of CPs: Residuals
        residual_X = Point3D_e1_PJCS[0] - Point3D_e2_PJCS[0]
        residual_Y = Point3D_e1_PJCS[1] - Point3D_e2_PJCS[1]
        residual_Z = Point3D_e1_PJCS[2] - Point3D_e2_PJCS[2]
        residual_P = np.linalg.norm(Point3D_e1_PJCS[0:2] - Point3D_e2_PJCS[0:2])
        residual_T = np.linalg.norm(Point3D_e1_PJCS[0:3] - Point3D_e2_PJCS[0:3])
        # analyse quality of CPs: Reprojection error
        RepErr_v2_e1 = MsTri.PointRepErrSigma(Point3D_e1_CLCS, Track1, cameras)
        RepErr_v3_e1 = MsTri.PointRepErrSigma(Point3D_e1_CLCS, Track2, cameras)
        RepErr_v2_e2 = MsTri.PointRepErrSigma(Point3D_e2_CLCS, Track2, cameras)
        RepErr_v3_e2 = MsTri.PointRepErrSigma(Point3D_e2_CLCS, Track1, cameras)
        # analyse quality of CPs: Triangulation angle
        TriAng_List = MsTri.PointTriAngList(Point3D_e1_CLCS, Track1, cameras) + \
                      MsTri.PointTriAngList(Point3D_e2_CLCS, Track2, cameras)
        TriAng_Max = max(TriAng_List)
        TriAng_Avg = np.mean(TriAng_List)
        TriAng_Min = min(TriAng_List)
        # analyse quality of CPs: epoch match number
        EpoMatch_Num = ViewNum1 * ViewNum2

        # construct result
        MarkerQua = [marker_id,
                     Point3D_e1_PJCS[0], Point3D_e1_PJCS[1], Point3D_e1_PJCS[2],
                     Point3D_e2_PJCS[0], Point3D_e2_PJCS[1], Point3D_e2_PJCS[2],
                     residual_X, residual_Y, residual_Z, residual_P, residual_T,
                     RepErr_v2_e1, RepErr_v3_e1, RepErr_v2_e2, RepErr_v3_e2,
                     TriAng_Max, TriAng_Avg, TriAng_Min, ViewNum1, ViewNum2, EpoMatch_Num]
        MarkersQua.append(MarkerQua)
    return MarkersQua


# 导出markers的质量结果
def exportMarkersQua(MarkersQua, path):
    '''
    input:
        MarkersQua = [MarkerQua,MarkerQua,...,MarkerQua]
            ,where MarkerQua = [marker_id,
                                    0
                                Point3D_e1[0], Point3D_e1[1], Point3D_e1[2], Point3D_e2[0], Point3D_e2[1], Point3D_e2[2],
                                      1              2              3             4              5               6
                                residual_X, residual_Y, residual_Z, residual_P, residual_T,
                                     7          8           9           10          11
                                RepErr_v2_e1, RepErr_v3_e1, RepErr_v2_e2, RepErr_v3_e2,
                                     12            13            14            15
                                TriAng_Max, TriAng_Avg, TriAng_Min, ViewNum1, ViewNum2, EpoMatch_Num]
                                     16         17          18         19        20          21
    output:
        _MarkersQua.txt
    '''
    reportFile = open(path, "w")
    fwriter = csv.writer(reportFile, delimiter='\t', lineterminator='\n')
    fwriter.writerow(['marker_id',
                      'X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2',
                      'ResX', 'ResY', 'ResZ', 'ResP', 'ResT',
                      'RepErr_v2_e1', 'RepErr_v3_e1', 'RepErr_v2_e2', 'RepErr_v3_e2',
                      'TriAng_Max', 'TriAng_Avg', 'TriAng_Min', 'ViewNum1', 'ViewNum2', 'EpoMatch_Num'])
    for MarkerQua in MarkersQua:
        fwriter.writerow(['{0}'.format(MarkerQua[0]),
                          '{0:0.5f}'.format(MarkerQua[1]),
                          '{0:0.5f}'.format(MarkerQua[2]),
                          '{0:0.5f}'.format(MarkerQua[3]),
                          '{0:0.5f}'.format(MarkerQua[4]),
                          '{0:0.5f}'.format(MarkerQua[5]),
                          '{0:0.5f}'.format(MarkerQua[6]),
                          '{0:0.5f}'.format(MarkerQua[7]),
                          '{0:0.5f}'.format(MarkerQua[8]),
                          '{0:0.5f}'.format(MarkerQua[9]),
                          '{0:0.5f}'.format(MarkerQua[10]),
                          '{0:0.5f}'.format(MarkerQua[11]),
                          '{0:0.3f}'.format(MarkerQua[12]),
                          '{0:0.3f}'.format(MarkerQua[13]),
                          '{0:0.3f}'.format(MarkerQua[14]),
                          '{0:0.3f}'.format(MarkerQua[15]),
                          '{0:0.3f}'.format(MarkerQua[16]),
                          '{0:0.3f}'.format(MarkerQua[17]),
                          '{0:0.3f}'.format(MarkerQua[18]),
                          '{0}'.format(MarkerQua[19]),
                          '{0}'.format(MarkerQua[20]),
                          '{0}'.format(MarkerQua[21])])
    reportFile.close()


# 对markers的质量进行统计并导出结果报告
def reportMarkersQua(MarkersQua, path):
    '''
    input:
        MarkersQua = [MarkerQua,MarkerQua,...,MarkerQua]
            ,where MarkerQua = [marker_id,
                                    0
                                Point3D_e1[0], Point3D_e1[1], Point3D_e1[2], Point3D_e2[0], Point3D_e2[1], Point3D_e2[2],
                                      1              2              3             4              5               6
                                residual_X, residual_Y, residual_Z, residual_P, residual_T,
                                     7          8           9           10          11
                                RepErr_v2_e1, RepErr_v3_e1, RepErr_v2_e2, RepErr_v3_e2,
                                     12            13            14            15
                                TriAng_Max, TriAng_Avg, TriAng_Min, ViewNum1, ViewNum2, EpoMatch_Num]
                                     16         17          18         19        20          21
    output:
        _MarkersQuaStatistic.txt
    '''
    # get statistic of residuals
    residual_X_List = [MarkerQua[7] for MarkerQua in MarkersQua]
    residual_Y_List = [MarkerQua[8] for MarkerQua in MarkersQua]
    residual_Z_List = [MarkerQua[9] for MarkerQua in MarkersQua]
    residual_P_List = [MarkerQua[10] for MarkerQua in MarkersQua]
    residual_T_List = [MarkerQua[11] for MarkerQua in MarkersQua]
    # RMSE
    residual_X_RMSE = FcSta.calculateRMSE(residual_X_List)
    residual_Y_RMSE = FcSta.calculateRMSE(residual_Y_List)
    residual_Z_RMSE = FcSta.calculateRMSE(residual_Z_List)
    residual_P_RMSE = FcSta.calculateRMSE(residual_P_List)
    residual_T_RMSE = FcSta.calculateRMSE(residual_T_List)
    # MAE
    residual_X_MAE = FcSta.calculateMAE(residual_X_List)
    residual_Y_MAE = FcSta.calculateMAE(residual_Y_List)
    residual_Z_MAE = FcSta.calculateMAE(residual_Z_List)
    residual_P_MAE = FcSta.calculateMAE(residual_P_List)
    residual_T_MAE = FcSta.calculateMAE(residual_T_List)
    # average
    residual_X_AVG = np.mean(residual_X_List)
    residual_Y_AVG = np.mean(residual_Y_List)
    residual_Z_AVG = np.mean(residual_Z_List)
    residual_P_AVG = np.mean(residual_P_List)
    residual_T_AVG = np.mean(residual_T_List)
    # standard deviation
    residual_X_STD = FcSta.calculateSTD(residual_X_List)
    residual_Y_STD = FcSta.calculateSTD(residual_Y_List)
    residual_Z_STD = FcSta.calculateSTD(residual_Z_List)
    residual_P_STD = FcSta.calculateSTD(residual_P_List)
    residual_T_STD = FcSta.calculateSTD(residual_T_List)

    # get distribution of quality item
    Distrib_residual_X = FcSta.listStatistic(residual_X_List)
    Distrib_residual_Y = FcSta.listStatistic(residual_Y_List)
    Distrib_residual_Z = FcSta.listStatistic(residual_Z_List)
    Distrib_residual_P = FcSta.listStatistic(residual_P_List)
    Distrib_residual_T = FcSta.listStatistic(residual_T_List)
    Distrib_RepErr_v2_e1 = FcSta.listStatistic([MarkerQua[12] for MarkerQua in MarkersQua])
    Distrib_RepErr_v3_e1 = FcSta.listStatistic([MarkerQua[13] for MarkerQua in MarkersQua])
    Distrib_RepErr_v2_e2 = FcSta.listStatistic([MarkerQua[14] for MarkerQua in MarkersQua])
    Distrib_RepErr_v3_e2 = FcSta.listStatistic([MarkerQua[15] for MarkerQua in MarkersQua])
    Distrib_TriAng_Max = FcSta.listStatistic([MarkerQua[16] for MarkerQua in MarkersQua])
    Distrib_TriAng_Avg = FcSta.listStatistic([MarkerQua[17] for MarkerQua in MarkersQua])
    Distrib_TriAng_Min = FcSta.listStatistic([MarkerQua[18] for MarkerQua in MarkersQua])
    Distrib_ViewNum1 = FcSta.listStatistic([MarkerQua[19] for MarkerQua in MarkersQua])
    Distrib_ViewNum2 = FcSta.listStatistic([MarkerQua[20] for MarkerQua in MarkersQua])
    Distrib_EpoMatch_Num = FcSta.listStatistic([MarkerQua[21] for MarkerQua in MarkersQua])

    # write file
    reportFile = open(path, "w")
    fwriter = csv.writer(reportFile, delimiter='\t', lineterminator='\n')
    # number of CTPs
    fwriter.writerow(['Res_X_RMSE(m)', 'Res_Y_RMSE(m)', 'Res_Z_RMSE(m)', 'Res_P_RMSE(m)', 'Res_T_RMSE(m)'])
    fwriter.writerow(['{0:0.5f}'.format(residual_X_RMSE),
                      '{0:0.5f}'.format(residual_Y_RMSE),
                      '{0:0.5f}'.format(residual_Z_RMSE),
                      '{0:0.5f}'.format(residual_P_RMSE),
                      '{0:0.5f}'.format(residual_T_RMSE)])
    fwriter.writerow(['Res_X_MAE(m)', 'Res_Y_MAE(m)', 'Res_Z_MAE(m)', 'Res_P_MAE(m)', 'Res_T_MAE(m)'])
    fwriter.writerow(['{0:0.5f}'.format(residual_X_MAE),
                      '{0:0.5f}'.format(residual_Y_MAE),
                      '{0:0.5f}'.format(residual_Z_MAE),
                      '{0:0.5f}'.format(residual_P_MAE),
                      '{0:0.5f}'.format(residual_T_MAE)])
    fwriter.writerow(['Res_X_AVG(m)', 'Res_Y_AVG(m)', 'Res_Z_AVG(m)', 'Res_P_AVG(m)', 'Res_T_AVG(m)'])
    fwriter.writerow(['{0:0.5f}'.format(residual_X_AVG),
                      '{0:0.5f}'.format(residual_Y_AVG),
                      '{0:0.5f}'.format(residual_Z_AVG),
                      '{0:0.5f}'.format(residual_P_AVG),
                      '{0:0.5f}'.format(residual_T_AVG)])
    fwriter.writerow(['Res_X_STD(m)', 'Res_Y_STD(m)', 'Res_Z_STD(m)', 'Res_P_STD(m)', 'Res_T_STD(m)'])
    fwriter.writerow(['{0:0.5f}'.format(residual_X_STD),
                      '{0:0.5f}'.format(residual_Y_STD),
                      '{0:0.5f}'.format(residual_Z_STD),
                      '{0:0.5f}'.format(residual_P_STD),
                      '{0:0.5f}'.format(residual_T_STD)])
    FcFile.writeListStatistic(fwriter, 'Distrib_residual_X', Distrib_residual_X)
    FcFile.writeListStatistic(fwriter, 'Distrib_residual_Y', Distrib_residual_Y)
    FcFile.writeListStatistic(fwriter, 'Distrib_residual_Z', Distrib_residual_Z)
    FcFile.writeListStatistic(fwriter, 'Distrib_residual_P', Distrib_residual_P)
    FcFile.writeListStatistic(fwriter, 'Distrib_residual_T', Distrib_residual_T)
    FcFile.writeListStatistic(fwriter, 'Distrib_RepErr_v2_e1', Distrib_RepErr_v2_e1)
    FcFile.writeListStatistic(fwriter, 'Distrib_RepErr_v2_e2', Distrib_RepErr_v2_e2)
    FcFile.writeListStatistic(fwriter, 'Distrib_RepErr_v3_e1', Distrib_RepErr_v3_e1)
    FcFile.writeListStatistic(fwriter, 'Distrib_RepErr_v3_e2', Distrib_RepErr_v3_e2)
    FcFile.writeListStatistic(fwriter, 'Distrib_TriAng_Max', Distrib_TriAng_Max)
    FcFile.writeListStatistic(fwriter, 'Distrib_TriAng_Avg', Distrib_TriAng_Avg)
    FcFile.writeListStatistic(fwriter, 'Distrib_TriAng_Min', Distrib_TriAng_Min)
    FcFile.writeListStatistic(fwriter, 'Distrib_ViewNum1', Distrib_ViewNum1)
    FcFile.writeListStatistic(fwriter, 'Distrib_ViewNum2', Distrib_ViewNum2)
    FcFile.writeListStatistic(fwriter, 'Distrib_EpoMatch_Num', Distrib_EpoMatch_Num)
    reportFile.close()
    return None


# 以GCP的【数据分析格式】导入为Metashape.Markers
def addMarkers_Analyse(chunk, camera_ids, Markers_Analyse, MarkerGroupName):
    '''
    get track, point and information from Metashape.Markers in document
    input:
        Markers_Analyse = [Marker,Marker,...,Marker]
            ,where Marker = [marker_id, marker_label, MarkerTrack]
            ,where MarkerTrack = [view,view,...,view]
            ,where view = [identifier, projection_info]
            ,where projection_info = [u, v]
    output:

    '''
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)
    cameras = chunk.cameras

    # creat new marker group
    MarkerGroup = chunk.addMarkerGroup()
    MarkerGroup.label = MarkerGroupName
    # creat new markers
    for Marker in Markers_Analyse:
        marker_id = Marker[0]
        marker_label = Marker[1]
        MarkerTrack = Marker[2]

        MarkerTrack_valid = []
        for view in MarkerTrack:
            # convert identifier to camera_id
            identifier = int(view[0])
            camera_id = camera_ids.get(identifier)
            if camera_id is not None:
                view[0] = camera_id
                MarkerTrack_valid.append(view)

        Point_Coord_CLCS = MsTri.Triangulate(MarkerTrack_valid, cameras)  # Vector([X,Y,Z])

        if Point_Coord_CLCS[-1] == 1:
            point_coord_CLCS = Metashape.Vector([Point_Coord_CLCS[0], Point_Coord_CLCS[1], Point_Coord_CLCS[2]])
            # point_coord_PJCS = MsTrans.transPointCoord_CLCS2PJCS(point_coord_CLCS, crs, M)
            # creat new marker
            marker = chunk.addMarker(point_coord_CLCS, visibility=False)
            # assign marker's attributions
            marker.group = MarkerGroup
            marker.label = 'Marker_{0}'.format(marker_label)
            # marker.reference.location = point_coord_PJCS
            marker.reference.accuracy = Metashape.Vector([10, 10, 10])  # in metres
            marker.reference.enabled = False
            # modify or append projection information, at the meanwhile, record the camera index list of MarkerTracks
            MarkerTracksCamera_IdList = []
            for view in MarkerTrack_valid:
                camera = cameras[view[0]]
                projection_info = view[1]
                projection_coord = Metashape.Vector([projection_info[0], projection_info[1]])
                # replace or append projection information of marker by MarkerTracks information
                marker.projections[camera] = Metashape.Marker.Projection(projection_coord, True)
                # record the camera index list of MarkerTracks
                MarkerTracksCamera_IdList.append(view[0])
            # delete projection of marker, which is not in MarkerTracks
            for camera, projection in marker.projections.items():
                if cameras.index(camera) not in MarkerTracksCamera_IdList:
                    del marker.projections[camera]


'''【数据分析格式】与【数据导入格式】相互转换'''

'''用于获取Markers的Metashape指标'''


# 导出Metashape对Markers计算的error指标
def getMarkerMsQua(chunk, MarkerGroupName, mode='both'):
    '''
    input:
        chunk
    output:
        MarkerErrors = [MarkerError,MarkerError,...,MarkerError]
            ,where MarkerError = [marker_id,marker_Enable,marker_label,
                                  Marker_Ref_X,Marker_Ref_Y,Marker_Ref_Z,
                                  Marker_Est_X,Marker_Est_Y,Marker_Est_Z,
                                  Marker_Acc_X,Marker_Acc_Y,Marker_Acc_Z,
                                  Marker_Err_X,Marker_Err_Y,Marker_Err_Z,
                                  error_meter, error_pixel,projectionsNum]
    '''
    MarkerErrors = []
    if MarkerGroupName:
        for markerGroup in chunk.marker_groups:
            if markerGroup.label != MarkerGroupName:
                continue
            # Now, we get the markerGroup of given MarkerGroupName
            for marker_id, marker in enumerate(chunk.markers):
                if marker.group != markerGroup:
                    continue
                if not MarkerMode(marker, mode):
                    continue
                MarkerError = calculateMarkerError(marker, chunk)
                MarkerError[0] = marker_id
                MarkerErrors.append(MarkerError)
    else:
        for marker_id, marker in enumerate(chunk.markers):
            if not MarkerMode(marker, mode):
                continue
            MarkerError = calculateMarkerError(marker, chunk)
            MarkerError[0] = marker_id
            MarkerErrors.append([MarkerError])
    return MarkerErrors


def exportMarkerErrors(MarkerErrors, path):
    '''
    input:
        MarkerErrors = [MarkerError,MarkerError,...,MarkerError]
            ,where MarkerError = [marker_id,marker_Enable,marker_label,
                                  Marker_Ref_X,Marker_Ref_Y,Marker_Ref_Z,
                                  Marker_Est_X,Marker_Est_Y,Marker_Est_Z,
                                  Marker_Acc_X,Marker_Acc_Y,Marker_Acc_Z,
                                  Marker_Err_X,Marker_Err_Y,Marker_Err_Z,
                                  error_meter, error_pixel,projectionsNum]
    output:
        .../MarkerErrors.txt
    '''
    # creat a file for writing results
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    for MarkerError in MarkerErrors:
        fwriter.writerow([MarkerError])
    File.close()


def importMarkerErrors(path):
    csv.field_size_limit(4000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        MarkerErrors = []
        for row in freader:
            MarkerErrors.append(eval(row[0]))
    return MarkerErrors


'''分析markers的质量'''


# 基于Metashape.Markers的【数据导入格式】分析markers的质量（包括Metashape的质量指标）
def getMarkersAllQua(chunk, Markers, camera_ids, CamerasEpoch, grid_id, MarkerFormat='Add', Triangulation='linear'):
    '''
    input:
        Markers = [Marker,Marker,...,Marker]    # MarkerFormat='Add'
            ,where Marker = [MarkerPoint,MarkerTrack,MarkerInform]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
            ,where MarkerTrack = [view,view,...,view]
                ,where view = [camera_id, projection_info]
                ,where projection_info = [u, v]
            ,where MarkerInform = [[point_id,Track_id],MarkerKey] or [[CrossMatch,CrossMatch,...,CrossMatch],MarkerKey]

        Markers = [Marker,Marker,...,Marker]    # MarkerFormat='Analyse'
            ,where Marker = [marker_id, marker_label, MarkerTrack]
            ,where MarkerTrack = [view,view,...,view]
            ,where view = [identifier, projection_info]
            ,where projection_info = [u, v]
        camera_ids = {identifier:camera_id}
        MarkerFormat='Add' or 'Analyse'
        Triangulation='nonlinear' or 'linear'
    output:
        MarkersQua = [MarkerQua,MarkerQua,...,MarkerQua]
            ,where MarkerQua = [id,
                                0
                                Point3D_e1_PJCS[0], Point3D_e1_PJCS[1], Point3D_e1_PJCS[2],
                                1                   2                   3
                                Point3D_e2_PJCS[0], Point3D_e2_PJCS[1], Point3D_e2_PJCS[2],
                                4                   5                   6
                                residual_X, residual_Y, residual_Z, residual_P, residual_T,
                                7           8           9           10          11
                                error_pixel, RepErr_v3_e1, RepErr_v3_e2, RepErr_v2_e1, RepErr_v2_e2,
                                12           13            14            15            16
                                TriAng_Max, TriAng_Avg, TriAng_Min,
                                17          18          19
                                ViewNum, ViewNum1, ViewNum2, EpoMatch_Num]
                                20       21        22        23
    '''
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)
    cameras = chunk.cameras

    MarkersAllQua = []
    progressbar = range(0, len(Markers), 10000)
    for i, Marker in enumerate(Markers):
        # get progress
        if i in progressbar and i != 0:
            print(i)

        # convert identifier to camera_id   $time cost: 2.75000000253082E-06 sec
        if MarkerFormat == 'Add':
            MarkerTrack = Marker[1]
        elif MarkerFormat == 'Analyse':
            MarkerTrack = Marker[2]
            for view in MarkerTrack:
                # convert identifier to camera_id
                camera_id = camera_ids[int(view[0])]
                view[0] = camera_id
        else:
            raise Exception("MarkerFormat should be 'Add' or 'Analyse'")

        # split MarkerTrack by epoch    $time cost: 3.28333333025436E-06 sec
        Track1 = []
        Track2 = []
        for view in MarkerTrack:
            # try to read the camera of this view, if the camera is not exist, this view will be skipped.
            camera_id = view[0]
            projection_info = view[1]
            if CamerasEpoch[camera_id] == 0:
                Track1.append([camera_id, projection_info])
            elif CamerasEpoch[camera_id] == 1:
                Track2.append([camera_id, projection_info])
            else:
                raise Exception("[Script]    the given epochs are not include all dates")

        # skip the case that sub-Track contains not enough images
        ViewNum = len(MarkerTrack)
        ViewNum1 = len(Track1)
        ViewNum2 = len(Track2)
        if ViewNum1 < 2 or ViewNum2 < 2:
            print('[Script]    Marker skipped:{0}'.format(i))
            continue

        # Triangulate   # $time cost: 0.000245050000001375 sec
        Point3D_e1_CLCS = MsTri.Triangulate(Track1, cameras)
        Point3D_e2_CLCS = MsTri.Triangulate(Track2, cameras)
        if Triangulation == 'nonlinear':
            # Triangulate with optimization   # $time cost: 0.0015783333333322 sec
            Point3D_e1_CLCS = MsTri.Triangulate_N(Point3D_e1_CLCS, Track1, cameras, WeightFunc=0)
            Point3D_e2_CLCS = MsTri.Triangulate_N(Point3D_e2_CLCS, Track2, cameras, WeightFunc=0)
            Point3D_e1_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e1_CLCS, crs, M)
            Point3D_e2_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e2_CLCS, crs, M)
        elif Triangulation == 'linear':
            Point3D_e1_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e1_CLCS, crs, M)
            Point3D_e2_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e2_CLCS, crs, M)
        else:
            raise Exception("Triangulation should be 'nonlinear' or 'linear'")

        # analyse quality of CPs: Residuals
        residual_X = Point3D_e1_PJCS[0] - Point3D_e2_PJCS[0]
        residual_Y = Point3D_e1_PJCS[1] - Point3D_e2_PJCS[1]
        residual_Z = Point3D_e1_PJCS[2] - Point3D_e2_PJCS[2]
        residual_P = np.linalg.norm(Point3D_e1_PJCS[0:2] - Point3D_e2_PJCS[0:2])
        residual_T = np.linalg.norm(Point3D_e1_PJCS[0:3] - Point3D_e2_PJCS[0:3])
        # analyse quality of CPs: Reprojection error   # $time cost: 0.0000469166666666372 sec
        RepErr_v2_e1 = MsTri.PointRepErrSigma(Point3D_e1_CLCS, Track1, cameras)
        RepErr_v2_e2 = MsTri.PointRepErrSigma(Point3D_e2_CLCS, Track2, cameras)
        # analyse quality of CPs: Reprojection error   # $time cost: 0.0000470333333358515 sec
        RepErr_v3_e1 = MsTri.PointRepErrSigma(Point3D_e1_CLCS, Track2, cameras)
        RepErr_v3_e2 = MsTri.PointRepErrSigma(Point3D_e2_CLCS, Track1, cameras)
        # analyse quality of CPs: Triangulation angle
        TriAng_List = MsTri.PointTriAngList(Point3D_e1_CLCS, Track1, cameras) + \
                      MsTri.PointTriAngList(Point3D_e2_CLCS, Track2, cameras)
        TriAng_Max = max(TriAng_List)
        TriAng_Avg = np.mean(TriAng_List)
        TriAng_Min = min(TriAng_List)
        # analyse quality of CPs: epoch match number
        EpoMatch_Num = ViewNum1 * ViewNum2
        # analyse quality of CPs: error pixel
        if Triangulation == 'nonlinear':
            Point_Coord_CLCS = MsTri.Triangulate(MarkerTrack, cameras)  # Vector([X,Y,Z])
            Point_Coord_CLCS = MsTri.Triangulate_N(Point_Coord_CLCS, MarkerTrack, cameras, WeightFunc=0)
            error_pixel = MsTri.PointRepErrSigma(Point_Coord_CLCS, MarkerTrack, cameras)
        elif Triangulation == 'linear':
            Point_Coord_CLCS = MsTri.Triangulate(MarkerTrack, cameras)  # Vector([X,Y,Z])
            error_pixel = MsTri.PointRepErrSigma(Point_Coord_CLCS, MarkerTrack, cameras)
        else:
            raise Exception("Triangulation should be 'nonlinear' or 'linear'")

        # use MarkerKey
        if grid_id:
            MarkerKey = getMarkerKey(grid_id, i)
        else:
            MarkerKey = i
        # construct result
        MarkerAllQua = [MarkerKey,
                        Point3D_e1_PJCS[0], Point3D_e1_PJCS[1], Point3D_e1_PJCS[2],
                        Point3D_e2_PJCS[0], Point3D_e2_PJCS[1], Point3D_e2_PJCS[2],
                        residual_X, residual_Y, residual_Z, residual_P, residual_T,
                        error_pixel, RepErr_v3_e1, RepErr_v3_e2, RepErr_v2_e1, RepErr_v2_e2,
                        TriAng_Max, TriAng_Avg, TriAng_Min,
                        ViewNum, ViewNum1, ViewNum2, EpoMatch_Num]
        MarkersAllQua.append(MarkerAllQua)
    return MarkersAllQua


# 导出markers的质量结果
def exportMarkersAllQua(MarkersAllQua, path):
    '''
    input:
        MarkersQua = [MarkerQua,MarkerQua,...,MarkerQua]
            ,where MarkerQua = [i,
                                0
                                Point3D_e1_PJCS[0], Point3D_e1_PJCS[1], Point3D_e1_PJCS[2],
                                1                   2                   3
                                Point3D_e2_PJCS[0], Point3D_e2_PJCS[1], Point3D_e2_PJCS[2],
                                4                   5                   6
                                residual_X, residual_Y, residual_Z, residual_P, residual_T,
                                7           8           9           10          11
                                error_pixel, RepErr_v3_e1, RepErr_v3_e2, RepErr_v2_e1, RepErr_v2_e2,
                                12           13            14            15            16
                                TriAng_Max, TriAng_Avg, TriAng_Min,
                                17          18          19
                                ViewNum, ViewNum1, ViewNum2, EpoMatch_Num]
                                20       21        22        23
    output:
        _MarkersAllQua.txt
    '''
    reportFile = open(path, "w")
    fwriter = csv.writer(reportFile, delimiter='\t', lineterminator='\n')
    fwriter.writerow(['i', 'X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2', 'ResX', 'ResY', 'ResZ', 'ResP', 'ResT',
                      'ErrorPixel', 'RepErr_v3_e1', 'RepErr_v3_e2', 'RepErr_v2_e1', 'RepErr_v2_e2',
                      'TriAng_Max', 'TriAng_Avg', 'TriAng_Min', 'ViewNum', 'ViewNum1', 'ViewNum2', 'EpoMatch_Num'])
    for MarkerAllQua in MarkersAllQua:
        fwriter.writerow(['{0}'.format(MarkerAllQua[0]),
                          '{0:0.5f}'.format(MarkerAllQua[1]),
                          '{0:0.5f}'.format(MarkerAllQua[2]),
                          '{0:0.5f}'.format(MarkerAllQua[3]),
                          '{0:0.5f}'.format(MarkerAllQua[4]),
                          '{0:0.5f}'.format(MarkerAllQua[5]),
                          '{0:0.5f}'.format(MarkerAllQua[6]),
                          '{0:0.5f}'.format(MarkerAllQua[7]),
                          '{0:0.5f}'.format(MarkerAllQua[8]),
                          '{0:0.5f}'.format(MarkerAllQua[9]),
                          '{0:0.5f}'.format(MarkerAllQua[10]),
                          '{0:0.5f}'.format(MarkerAllQua[11]),
                          '{0:0.5f}'.format(MarkerAllQua[12]),
                          '{0:0.5f}'.format(MarkerAllQua[13]),
                          '{0:0.5f}'.format(MarkerAllQua[14]),
                          '{0:0.5f}'.format(MarkerAllQua[15]),
                          '{0:0.5f}'.format(MarkerAllQua[16]),
                          '{0:0.3f}'.format(MarkerAllQua[17]),
                          '{0:0.3f}'.format(MarkerAllQua[18]),
                          '{0:0.3f}'.format(MarkerAllQua[19]),
                          '{0}'.format(MarkerAllQua[20]),
                          '{0}'.format(MarkerAllQua[21]),
                          '{0}'.format(MarkerAllQua[22]),
                          '{0}'.format(MarkerAllQua[23])])
    reportFile.close()


# 对markers的质量进行统计并导出结果报告
def reportMarkersAllQua(MarkersAllQua, path):
    '''
    input:
        MarkersQua = [MarkerQua,MarkerQua,...,MarkerQua]
            ,where MarkerQua = [i,
                                0
                                Point3D_e1_PJCS[0], Point3D_e1_PJCS[1], Point3D_e1_PJCS[2],
                                1                   2                   3
                                Point3D_e2_PJCS[0], Point3D_e2_PJCS[1], Point3D_e2_PJCS[2],
                                4                   5                   6
                                residual_X, residual_Y, residual_Z, residual_P, residual_T,
                                7           8           9           10          11
                                error_pixel, RepErr_v3_e1, RepErr_v3_e2, RepErr_v2_e1, RepErr_v2_e2,
                                12           13            14            15            16
                                TriAng_Max, TriAng_Avg, TriAng_Min,
                                17          18          19
                                ViewNum, ViewNum1, ViewNum2, EpoMatch_Num]
                                20       21        22        23
    output:
        _MarkersQuaSta.txt
    '''
    # get statistic of residuals
    residual_X_List = [AllQua[7] for AllQua in MarkersAllQua]
    residual_Y_List = [AllQua[8] for AllQua in MarkersAllQua]
    residual_Z_List = [AllQua[9] for AllQua in MarkersAllQua]
    residual_P_List = [AllQua[10] for AllQua in MarkersAllQua]
    residual_T_List = [AllQua[11] for AllQua in MarkersAllQua]
    # RMSE
    residual_X_RMSE = FcSta.calculateRMSE(residual_X_List)
    residual_Y_RMSE = FcSta.calculateRMSE(residual_Y_List)
    residual_Z_RMSE = FcSta.calculateRMSE(residual_Z_List)
    residual_P_RMSE = FcSta.calculateRMSE(residual_P_List)
    residual_T_RMSE = FcSta.calculateRMSE(residual_T_List)
    # MAE
    residual_X_MAE = FcSta.calculateMAE(residual_X_List)
    residual_Y_MAE = FcSta.calculateMAE(residual_Y_List)
    residual_Z_MAE = FcSta.calculateMAE(residual_Z_List)
    residual_P_MAE = FcSta.calculateMAE(residual_P_List)
    residual_T_MAE = FcSta.calculateMAE(residual_T_List)
    # average
    residual_X_AVG = np.mean(residual_X_List)
    residual_Y_AVG = np.mean(residual_Y_List)
    residual_Z_AVG = np.mean(residual_Z_List)
    residual_P_AVG = np.mean(residual_P_List)
    residual_T_AVG = np.mean(residual_T_List)
    # standard deviation
    residual_X_STD = FcSta.calculateSTD(residual_X_List)
    residual_Y_STD = FcSta.calculateSTD(residual_Y_List)
    residual_Z_STD = FcSta.calculateSTD(residual_Z_List)
    residual_P_STD = FcSta.calculateSTD(residual_P_List)
    residual_T_STD = FcSta.calculateSTD(residual_T_List)

    # get distribution of quality item
    Distrib_residual_X = FcSta.listStatistic(residual_X_List)
    Distrib_residual_Y = FcSta.listStatistic(residual_Y_List)
    Distrib_residual_Z = FcSta.listStatistic(residual_Z_List)
    Distrib_residual_P = FcSta.listStatistic(residual_P_List)
    Distrib_residual_T = FcSta.listStatistic(residual_T_List)
    Distrib_ErrorPixel = FcSta.listStatistic([AllQua[12] for AllQua in MarkersAllQua])
    Distrib_RepErr_v3_e1 = FcSta.listStatistic([AllQua[13] for AllQua in MarkersAllQua])
    Distrib_RepErr_v3_e2 = FcSta.listStatistic([AllQua[14] for AllQua in MarkersAllQua])
    Distrib_RepErr_v2_e1 = FcSta.listStatistic([AllQua[15] for AllQua in MarkersAllQua])
    Distrib_RepErr_v2_e2 = FcSta.listStatistic([AllQua[16] for AllQua in MarkersAllQua])
    Distrib_TriAng_Max = FcSta.listStatistic([AllQua[17] for AllQua in MarkersAllQua])
    Distrib_TriAng_Avg = FcSta.listStatistic([AllQua[18] for AllQua in MarkersAllQua])
    Distrib_TriAng_Min = FcSta.listStatistic([AllQua[19] for AllQua in MarkersAllQua])
    Distrib_ViewNum = FcSta.listStatistic([AllQua[20] for AllQua in MarkersAllQua])
    Distrib_ViewNum1 = FcSta.listStatistic([AllQua[21] for AllQua in MarkersAllQua])
    Distrib_ViewNum2 = FcSta.listStatistic([AllQua[22] for AllQua in MarkersAllQua])
    Distrib_EpoMatch_Num = FcSta.listStatistic([AllQua[23] for AllQua in MarkersAllQua])

    # write file
    reportFile = open(path, "w")
    fwriter = csv.writer(reportFile, delimiter='\t', lineterminator='\n')
    # number of CTPs
    fwriter.writerow(['Res_X_RMSE(m)', 'Res_Y_RMSE(m)', 'Res_Z_RMSE(m)', 'Res_P_RMSE(m)', 'Res_T_RMSE(m)'])
    fwriter.writerow(['{0:0.5f}'.format(residual_X_RMSE),
                      '{0:0.5f}'.format(residual_Y_RMSE),
                      '{0:0.5f}'.format(residual_Z_RMSE),
                      '{0:0.5f}'.format(residual_P_RMSE),
                      '{0:0.5f}'.format(residual_T_RMSE)])
    fwriter.writerow(['Res_X_MAE(m)', 'Res_Y_MAE(m)', 'Res_Z_MAE(m)', 'Res_P_MAE(m)', 'Res_T_MAE(m)'])
    fwriter.writerow(['{0:0.5f}'.format(residual_X_MAE),
                      '{0:0.5f}'.format(residual_Y_MAE),
                      '{0:0.5f}'.format(residual_Z_MAE),
                      '{0:0.5f}'.format(residual_P_MAE),
                      '{0:0.5f}'.format(residual_T_MAE)])
    fwriter.writerow(['Res_X_AVG(m)', 'Res_Y_AVG(m)', 'Res_Z_AVG(m)', 'Res_P_AVG(m)', 'Res_T_AVG(m)'])
    fwriter.writerow(['{0:0.5f}'.format(residual_X_AVG),
                      '{0:0.5f}'.format(residual_Y_AVG),
                      '{0:0.5f}'.format(residual_Z_AVG),
                      '{0:0.5f}'.format(residual_P_AVG),
                      '{0:0.5f}'.format(residual_T_AVG)])
    fwriter.writerow(['Res_X_STD(m)', 'Res_Y_STD(m)', 'Res_Z_STD(m)', 'Res_P_STD(m)', 'Res_T_STD(m)'])
    fwriter.writerow(['{0:0.5f}'.format(residual_X_STD),
                      '{0:0.5f}'.format(residual_Y_STD),
                      '{0:0.5f}'.format(residual_Z_STD),
                      '{0:0.5f}'.format(residual_P_STD),
                      '{0:0.5f}'.format(residual_T_STD)])
    FcFile.writeListStatistic(fwriter, 'Distrib_residual_X', Distrib_residual_X)
    FcFile.writeListStatistic(fwriter, 'Distrib_residual_Y', Distrib_residual_Y)
    FcFile.writeListStatistic(fwriter, 'Distrib_residual_Z', Distrib_residual_Z)
    FcFile.writeListStatistic(fwriter, 'Distrib_residual_P', Distrib_residual_P)
    FcFile.writeListStatistic(fwriter, 'Distrib_residual_T', Distrib_residual_T)
    FcFile.writeListStatistic(fwriter, 'Distrib_ErrorPixel', Distrib_ErrorPixel)
    FcFile.writeListStatistic(fwriter, 'Distrib_RepErr_v3_e1', Distrib_RepErr_v3_e1)
    FcFile.writeListStatistic(fwriter, 'Distrib_RepErr_v3_e2', Distrib_RepErr_v3_e2)
    FcFile.writeListStatistic(fwriter, 'Distrib_RepErr_v2_e1', Distrib_RepErr_v2_e1)
    FcFile.writeListStatistic(fwriter, 'Distrib_RepErr_v2_e2', Distrib_RepErr_v2_e2)
    FcFile.writeListStatistic(fwriter, 'Distrib_TriAng_Max', Distrib_TriAng_Max)
    FcFile.writeListStatistic(fwriter, 'Distrib_TriAng_Avg', Distrib_TriAng_Avg)
    FcFile.writeListStatistic(fwriter, 'Distrib_TriAng_Min', Distrib_TriAng_Min)
    FcFile.writeListStatistic(fwriter, 'Distrib_ViewNum', Distrib_ViewNum)
    FcFile.writeListStatistic(fwriter, 'Distrib_ViewNum1', Distrib_ViewNum1)
    FcFile.writeListStatistic(fwriter, 'Distrib_ViewNum2', Distrib_ViewNum2)
    FcFile.writeListStatistic(fwriter, 'Distrib_EpoMatch_Num', Distrib_EpoMatch_Num)
    reportFile.close()
    return None


def readMarkersAllQua(path):
    '''
    input:
        Qua.txt
    output:
        MarkersQua = [Qua,Qua,...,Qua]
            Qua = [i,
                   0
                   Point3D_e1_PJCS[0], Point3D_e1_PJCS[1], Point3D_e1_PJCS[2],
                   1                   2                   3
                   Point3D_e2_PJCS[0], Point3D_e2_PJCS[1], Point3D_e2_PJCS[2],
                   4                   5                   6
                   residual_X, residual_Y, residual_Z, residual_P, residual_T,
                   7           8           9           10          11
                   error_pixel, RepErr_v3_e1, RepErr_v3_e2, RepErr_v2_e1, RepErr_v2_e2,
                   12           13            14            15            16
                   TriAng_Max, TriAng_Avg, TriAng_Min,
                   17          18          19
                   ViewNum, ViewNum1, ViewNum2, EpoMatch_Num]
                   20       21        22        23
    '''
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        head = next(freader)  # skip header
        MarkersQua = []
        for row in freader:
            MarkersQua.append([eval(item) for item in row])
    return MarkersQua


def readMarkersAllQuaReport(path):
    '''
    input:
        QuaSta.txt
    output:
        QuaSta
    Res_X_RMSE(m)   Res_Y_RMSE(m)   Res_Z_RMSE(m)   Res_P_RMSE(m)   Res_T_RMSE(m)
        0,1             2,3             4,5             6,7             8,9
    Res_X_MAE(m)    Res_Y_MAE(m)    Res_Z_MAE(m)    Res_P_MAE(m)    Res_T_MAE(m)
       10,11           12,13           14,15           16,17           18,19
    Res_X_AVG(m)    Res_Y_AVG(m)    Res_Z_AVG(m)    Res_P_AVG(m)    Res_T_AVG(m)
       20,21           22,23           24,25           26,27           28,29
    Res_X_STD(m)    Res_Y_STD(m)    Res_Z_STD(m)    Res_P_STD(m)    Res_T_STD(m)
       30,31           32,33           34,35           36,37           38,39
    Distrib_residual_X_AVG      Distrib_residual_X_MIN      Distrib_residual_X_25P      Distrib_residual_X_50P      Distrib_residual_X_75P      Distrib_residual_X_MAX
    40,41                       42,43                       44,45                       46,47                       48,49                       50,51
    Distrib_residual_Z_AVG      Distrib_residual_Z_MIN      Distrib_residual_Z_25P      Distrib_residual_Z_50P      Distrib_residual_Z_75P      Distrib_residual_Z_MAX
    52,53
    Distrib_residual_Y_AVG      Distrib_residual_Y_MIN      Distrib_residual_Y_25P      Distrib_residual_Y_50P      Distrib_residual_Y_75P      Distrib_residual_Y_MAX
    64,65
    Distrib_residual_P_AVG      Distrib_residual_P_MIN      Distrib_residual_P_25P      Distrib_residual_P_50P      Distrib_residual_P_75P      Distrib_residual_P_MAX
    76,77
    Distrib_residual_T_AVG      Distrib_residual_T_MIN      Distrib_residual_T_25P      Distrib_residual_T_50P      Distrib_residual_T_75P      Distrib_residual_T_MAX
    88,89
    Distrib_error_pixel_AVG     Distrib_error_pixel_MIN     Distrib_error_pixel_25P     Distrib_error_pixel_50P     Distrib_error_pixel_75P     Distrib_error_pixel_MAX
    100,101
    Distrib_RepErr_v3_e1_AVG    Distrib_RepErr_v3_e1_MIN    Distrib_RepErr_v3_e1_25P    Distrib_RepErr_v3_e1_50P    Distrib_RepErr_v3_e1_75P    Distrib_RepErr_v3_e1_MAX
    112,113
    Distrib_RepErr_v3_e2_AVG    Distrib_RepErr_v3_e2_MIN    Distrib_RepErr_v3_e2_25P    Distrib_RepErr_v3_e2_50P    Distrib_RepErr_v3_e2_75P    Distrib_RepErr_v3_e2_MAX
    124,125
    Distrib_RepErr_v2_e1_AVG    Distrib_RepErr_v2_e1_MIN    Distrib_RepErr_v2_e1_25P    Distrib_RepErr_v2_e1_50P    Distrib_RepErr_v2_e1_75P    Distrib_RepErr_v2_e1_MAX
    136,137
    Distrib_RepErr_v2_e2_AVG    Distrib_RepErr_v2_e2_MIN    Distrib_RepErr_v2_e2_25P    Distrib_RepErr_v2_e2_50P    Distrib_RepErr_v2_e2_75P    Distrib_RepErr_v2_e2_MAX
    148,149
    Distrib_TriAng_Max_AVG      Distrib_TriAng_Max_MIN      Distrib_TriAng_Max_25P      Distrib_TriAng_Max_50P      Distrib_TriAng_Max_75P      Distrib_TriAng_Max_MAX
    160,161
    Distrib_TriAng_Avg_AVG      Distrib_TriAng_Avg_MIN      Distrib_TriAng_Avg_25P      Distrib_TriAng_Avg_50P      Distrib_TriAng_Avg_75P      Distrib_TriAng_Avg_MAX
    172,173
    Distrib_TriAng_Min_AVG      Distrib_TriAng_Min_MIN      Distrib_TriAng_Min_25P      Distrib_TriAng_Min_50P      Distrib_TriAng_Min_75P      Distrib_TriAng_Min_MAX
    184,185
    Distrib_ViewNum_AVG         Distrib_ViewNum_MIN         Distrib_ViewNum_25P         Distrib_ViewNum_50P         Distrib_ViewNum_75P         Distrib_ViewNum_MAX
    196,197
    Distrib_ViewNum1_AVG        Distrib_ViewNum1_MIN        Distrib_ViewNum1_25P        Distrib_ViewNum1_50P        Distrib_ViewNum1_75P        Distrib_ViewNum1_MAX
    208,209
    Distrib_ViewNum2_AVG        Distrib_ViewNum2_MIN        Distrib_ViewNum2_25P        Distrib_ViewNum2_50P        Distrib_ViewNum2_75P        Distrib_ViewNum2_MAX
    220,221
    Distrib_EpoMatch_Num_AVG    Distrib_EpoMatch_Num_MIN    Distrib_EpoMatch_Num_25P    Distrib_EpoMatch_Num_50P    Distrib_EpoMatch_Num_75P    Distrib_EpoMatch_Num_MAX
    232,233
    '''
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        fileData = []
        MarkersQuaSta = {}
        for row in freader:
            fileData.append(row)
        for i in range(0, 8, 2):
            for j in range(len(fileData[i])):
                MarkersQuaSta[fileData[i][j]] = eval(fileData[i + 1][j])
        for i in range(8, len(fileData), 2):
            MarkersQuaSta[fileData[i][0][0:-4] + '_AVG'] = eval(fileData[i + 1][0])
            MarkersQuaSta[fileData[i][0][0:-4] + '_MIN'] = eval(fileData[i + 1][1])
            MarkersQuaSta[fileData[i][0][0:-4] + '_25P'] = eval(fileData[i + 1][2])
            MarkersQuaSta[fileData[i][0][0:-4] + '_50P'] = eval(fileData[i + 1][3])
            MarkersQuaSta[fileData[i][0][0:-4] + '_75P'] = eval(fileData[i + 1][4])
            MarkersQuaSta[fileData[i][0][0:-4] + '_MAX'] = eval(fileData[i + 1][5])
    return MarkersQuaSta


'''只分析markers的error pixel'''


def getMarkersRes(chunk, Markers, camera_ids, CamerasEpoch,
                  Index='ErrorPixel', MarkerFormat='Add', Triangulation='linear'):
    '''
    这个函数只适用于SelectCTPs_Optimal_pairwise_part2.py
    input:
        Index = 'ErrorPixel' or 'Residuals' or 'ReprojError' or 'TriAngles'
        Markers = [Marker,Marker,...,Marker]    # MarkerFormat='Add'
            ,where Marker = [MarkerPoint,MarkerTrack,MarkerInform]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
            ,where MarkerTrack = [view,view,...,view]
                ,where view = [camera_id, projection_info]
                ,where projection_info = [u, v]
            ,where MarkerInform = [[point_id,Track_id],MarkerKey] or [[CrossMatch,CrossMatch,...,CrossMatch],MarkerKey]

        Markers = [Marker,Marker,...,Marker]    # MarkerFormat='Analyse'
            ,where Marker = [marker_id, marker_label, MarkerTrack]
            ,where MarkerTrack = [view,view,...,view]
            ,where view = [identifier, projection_info]
            ,where projection_info = [u, v]
        camera_ids = {identifier:camera_id}
        MarkerFormat='Add' or 'Analyse'
        Triangulation='nonlinear' or 'linear'
    output:
        MarkersRes = [MarkerRes,MarkerRes,...,MarkerRes]
            ,if Index == 'ErrorPixel': MarkerRes = [MarkerKey, error_pixel, ViewNum, ViewNum1, ViewNum2]
            ,if Index == 'Residuals': MarkerRes = [MarkerKey, RX, RY, RZ, RP, RT, ViewNum, ViewNum1, ViewNum2]
            ,if Index == 'ReprojError': MarkerRes = [MarkerKey, ReprojError, ViewNum, ViewNum1, ViewNum2]
            ,if Index == 'TriAngles': MarkerRes = [MarkerKey, TriAng_Avg, ViewNum, ViewNum1, ViewNum2]
            ,if Index == 'MixReprojError': MarkerRes = [MarkerKey, MixRepErr, ViewNum, ViewNum1, ViewNum2]
    '''
    cameras = chunk.cameras
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)

    MarkersRes = []
    for i, Marker in enumerate(Markers):
        # combine marker grid_id and index in grid for packed to numpy ndarray
        MarkerKey = Marker[2][-1]

        # convert identifier to camera_id   $time cost: 2.75000000253082E-06 sec
        if MarkerFormat == 'Add':
            MarkerTrack = Marker[1]
        elif MarkerFormat == 'Analyse':
            MarkerTrack = Marker[2]
            for view in MarkerTrack:
                # convert identifier to camera_id
                camera_id = camera_ids[int(view[0])]
                view[0] = camera_id
        else:
            raise Exception("MarkerFormat should be 'Add' or 'Analyse'")

        # split MarkerTrack by epoch    $time cost: 3.28333333025436E-06 sec
        Track1 = []
        Track2 = []
        for view in MarkerTrack:
            # try to read the camera of this view, if the camera is not exist, this view will be skipped.
            camera_id = view[0]
            projection_info = view[1]
            if CamerasEpoch[camera_id] == 0:
                Track1.append([camera_id, projection_info])
            elif CamerasEpoch[camera_id] == 1:
                Track2.append([camera_id, projection_info])
            else:
                raise Exception("[Script]    the given epochs are not include all dates")
        ViewNum = len(MarkerTrack)
        ViewNum1 = len(Track1)
        ViewNum2 = len(Track2)

        # Calculates the specified metric
        if Index == 'ErrorPixel':
            # analyse quality of CPs: error pixel
            if Triangulation == 'nonlinear':
                Point_Coord_CLCS = MsTri.Triangulate(MarkerTrack, cameras)  # Vector([X,Y,Z])
                Point_Coord_CLCS = MsTri.Triangulate_N(Point_Coord_CLCS, MarkerTrack, cameras, WeightFunc=0)
                error_pixel = MsTri.PointRepErrSigma(Point_Coord_CLCS, MarkerTrack, cameras)  # error_pixel
            elif Triangulation == 'linear':
                Point_Coord_CLCS = MsTri.Triangulate(MarkerTrack, cameras)  # Vector([X,Y,Z])
                error_pixel = MsTri.PointRepErrSigma(Point_Coord_CLCS, MarkerTrack, cameras)
            else:
                raise Exception("Triangulation should be 'nonlinear' or 'linear'")
            # construct result
            MarkersRes.append([MarkerKey, error_pixel, ViewNum, ViewNum1, ViewNum2])
        else:
            # skip the case that sub-Track contains not enough images
            if ViewNum1 < 2 or ViewNum2 < 2:
                print('[Script]    Marker skipped')
                continue
            # Triangulate   # $time cost: 0.000245050000001375 sec
            Point3D_e1_CLCS = MsTri.Triangulate(Track1, cameras)
            Point3D_e2_CLCS = MsTri.Triangulate(Track2, cameras)
            if Index == 'Residuals':
                if Triangulation == 'nonlinear':
                    # Triangulate with optimization   # $time cost: 0.0015783333333322 sec
                    Point3D_e1_CLCS = MsTri.Triangulate_N(Point3D_e1_CLCS, Track1, cameras, WeightFunc=0)
                    Point3D_e2_CLCS = MsTri.Triangulate_N(Point3D_e2_CLCS, Track2, cameras, WeightFunc=0)
                    Point3D_e1_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e1_CLCS, crs, M)
                    Point3D_e2_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e2_CLCS, crs, M)
                elif Triangulation == 'linear':
                    Point3D_e1_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e1_CLCS, crs, M)
                    Point3D_e2_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e2_CLCS, crs, M)
                else:
                    raise Exception("Triangulation should be 'nonlinear' or 'linear'")
                RX = Point3D_e1_PJCS[0] - Point3D_e2_PJCS[0]
                RY = Point3D_e1_PJCS[1] - Point3D_e2_PJCS[1]
                RZ = Point3D_e1_PJCS[2] - Point3D_e2_PJCS[2]
                RP = np.linalg.norm(Point3D_e1_PJCS[0:2] - Point3D_e2_PJCS[0:2])
                RT = np.linalg.norm(Point3D_e1_PJCS[0:3] - Point3D_e2_PJCS[0:3])
                # construct result
                MarkersRes.append([MarkerKey, RX, RY, RZ, RP, RT, ViewNum, ViewNum1, ViewNum2])
            elif Index == 'ReprojError':
                # analyse quality of CPs: Reprojection error   # $time cost: 0.0000470333333358515 sec
                RepErr_v3_e1 = MsTri.PointRepErrSigma(Point3D_e1_CLCS, Track2, cameras)
                RepErr_v3_e2 = MsTri.PointRepErrSigma(Point3D_e2_CLCS, Track1, cameras)
                # construct result
                MarkersRes.append([MarkerKey, (RepErr_v3_e1 + RepErr_v3_e2) / 2, ViewNum, ViewNum1, ViewNum2])
            elif Index == 'MixReprojError':
                # analyse quality of CPs: error pixel
                if Triangulation == 'nonlinear':
                    Point_Coord_CLCS = MsTri.Triangulate(MarkerTrack, cameras)  # Vector([X,Y,Z])
                    Point_Coord_CLCS = MsTri.Triangulate_N(Point_Coord_CLCS, MarkerTrack, cameras, WeightFunc=0)
                    error_pixel = MsTri.PointRepErrSigma(Point_Coord_CLCS, MarkerTrack, cameras)  # error_pixel
                elif Triangulation == 'linear':
                    Point_Coord_CLCS = MsTri.Triangulate(MarkerTrack, cameras)  # Vector([X,Y,Z])
                    error_pixel = MsTri.PointRepErrSigma(Point_Coord_CLCS, MarkerTrack, cameras)
                else:
                    raise Exception("Triangulation should be 'nonlinear' or 'linear'")
                # analyse quality of CPs: Reprojection error   # $time cost: 0.0000470333333358515 sec
                RepErr_v3_e1 = MsTri.PointRepErrSigma(Point3D_e1_CLCS, Track2, cameras)
                RepErr_v3_e2 = MsTri.PointRepErrSigma(Point3D_e2_CLCS, Track1, cameras)
                # construct result
                MarkersRes.append(
                    [MarkerKey, (error_pixel * 2 + RepErr_v3_e1 + RepErr_v3_e2) / 4, ViewNum, ViewNum1, ViewNum2])
            elif Index == 'TriAngles':
                # analyse quality of CPs: Triangulation angle
                TriAng_List = MsTri.PointTriAngList(Point3D_e1_CLCS, Track1, cameras) + \
                              MsTri.PointTriAngList(Point3D_e2_CLCS, Track2, cameras)
                TriAng_Avg = np.mean(TriAng_List)
                # construct result
                MarkersRes.append([MarkerKey, TriAng_Avg, ViewNum, ViewNum1, ViewNum2])
    return MarkersRes


def exportMarkerResiduals(path, MarkerResiduals):
    reportFile = open(path, "w")
    fwriter = csv.writer(reportFile, delimiter='\t', lineterminator='\n')
    fwriter.writerow(['Error(pix)', 'RepErr_v3_e1', 'RepErr_v3_e2'])
    for MksRes in MarkerResiduals:
        fwriter.writerow(['{0:0.4f}'.format(MksRes[0]), '{0:0.4f}'.format(MksRes[1]), '{0:0.4f}'.format(MksRes[2])])
    reportFile.close()


########################################################################################################################
def getMarkerTrack(marker, cameras, Index='identifier'):
    '''
    get track from Metashape.Marker
    '''
    MarkerTrack = []
    for camera, projection in marker.projections.items():
        if camera.type == Metashape.Camera.Type.Keyframe:
            print('[Script]    the sensor is keyframe!:', camera.label)
            continue  # skipping Keyframes
        elif not camera.transform:
            print('[Script]    camera has no transform!:', camera.label)
            continue  # skipping no transform
        elif not camera.enabled:
            print('[Script]    camera is not enabled!', camera.label)
            continue  # skipping camera not used
        elif not projection.valid:  # 暂时描述projection是否被删除，直到被确认后才会被删除
            continue
        elif not projection.pinned:  # False:灰色小旗的点，True: 蓝色和绿色
            continue
        if Index == 'identifier':
            identifier = ConnectData_Metashape.getCameraIdentifier(camera)
            view = [identifier, [projection.coord[0], projection.coord[1]]]
        elif Index == 'camera_id':
            camera_id = cameras.index(camera)
            view = [camera_id, [projection.coord[0], projection.coord[1]]]
        MarkerTrack.append(view)
    return MarkerTrack


def MarkerMode(marker, mode):
    '''
    judge whether the marker conform to the given mode
    '''
    if mode == 'check':
        if marker.reference.enabled:
            return 1
        else:
            return 0
    elif mode == 'uncheck':
        if not marker.reference.enabled:
            return 1
        else:
            return 0
    elif mode == 'both':
        return 1
    else:
        print('wrong mode!')
        return 0


def calculateMarkerError(marker, chunk):
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)
    cameras = chunk.cameras

    Marker_Ref = marker.reference.location  # Vector([X,Y,Z])
    MarkerTrack = getMarkerTrack(marker, cameras, Index='camera_id')
    Point_Coord_CLCS = MsTri.Triangulate(MarkerTrack, cameras)  # Vector([X,Y,Z])
    Point_Coord_CLCS = MsTri.Triangulate_N(Point_Coord_CLCS, MarkerTrack, cameras, WeightFunc=0)
    Marker_Est = MsTrans.transPointCoord_CLCS2PJCS(Point_Coord_CLCS, crs, M)  # Vector([X,Y,Z])
    Marker_Acc = marker.reference.accuracy  # Vector([X,Y,Z])
    Marker_Err = Marker_Est - Marker_Ref  # Vector([X,Y,Z])
    error_meter = Marker_Err.norm()
    error_pixel = MsTri.PointRepErrSigma(Point_Coord_CLCS, MarkerTrack, cameras)
    return [-1, marker.reference.enabled, marker.label,
            Marker_Ref[0], Marker_Ref[1], Marker_Ref[2],
            Marker_Est[0], Marker_Est[1], Marker_Est[2],
            Marker_Acc[0], Marker_Acc[1], Marker_Acc[2],
            Marker_Err[0], Marker_Err[1], Marker_Err[2],
            error_meter, error_pixel, len(MarkerTrack)]


def getMarkerKey(grid_id, i):
    # a,b,i不超过6位数
    a = grid_id[0]
    b = grid_id[1]
    return i * 1e+12 + a * 1e+6 + b


def readMarkerKey(MarkerKey):
    # a,b,i不超过6位数
    (i, ab) = divmod(MarkerKey, 1e+12)
    grid_id = divmod(ab, 1e+6)
    return [grid_id, i]
