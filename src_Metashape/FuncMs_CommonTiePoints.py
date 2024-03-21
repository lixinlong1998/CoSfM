import math
import csv
import numpy as np
from scipy.stats import gaussian_kde
from src_CFTM import Func_SpatialDistribution
from src_CFTM import Func_CommonTiePoints
from src_CFTM import Func_Triangulate as Tri
from src_CFTM import Func_Statistic as FcSta
from src_CFTM import Func_Files as FcFile
import src_Metashape.FuncMs_Transform as MsTrans


# def getICTPs(chunk, Points, ICTPs_IdList, ICTPsTracks, CTPsTyp, CTPsDen, CTPsQua, CTPsRetriCoord):
#     '''
#     input:
#         chunk =
#         Points={point_id:Point,...}
#             ,where Point = [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, track_id]
#                    Coordinates=[[X, Y, Z],[X, Y, Z]]   # in CLCS & PJCS respectively
#                    Covariance = 3x3 covariance matrix
#                    StandardDeviation = [stdX, stdY, stdZ]  # in mm
#                    Quality = [RE, RU, PA, IC]
#         ICTPs_IdList =
#         ICTPsTracks =
#         CTPsTyp = [ETP_Type,ETP_Type,...,ETP_Type]
#         CTPsDen = [ETPDensity, ETPDensity, ..., ETPDensity]
#             ,where ETPDensity = [point_id, density]
#         CTPsQua = [CTPQua,CTPQua,...,CTPQua]
#             ,where CTPQua = [point_id,
#                              0
#                              Point3D_e1_PJCS[0], Point3D_e1_PJCS[1], Point3D_e1_PJCS[2],
#                              1                   2                   3
#                              Point3D_e2_PJCS[0], Point3D_e2_PJCS[1], Point3D_e2_PJCS[2],
#                              4                   5                   6
#                              residual_X, residual_Y, residual_Z, residual_P, residual_T,
#                              7           8           9           10          11
#                              error_pixel, RepErr_v3_e1, RepErr_v3_e2, RepErr_v2_e1, RepErr_v2_e2,
#                              12           13            14            15            16
#                              TriAng_Max, TriAng_Avg, TriAng_Min,
#                              17          18          19
#                              ViewNum, ViewNum1, ViewNum2, EpoMatch_Num]
#                              20       21        22        23
#         CTPsRetriCoord = [Point3D,Point3D,...,Point3D]
#             ,where Point3D = ndarray([X,Y,Z,1])
#     output:
#         ICTPsData = [MCTPData,MCTPData,...,MCTPData]
#             ,where MCTPData = [point_id, track_id, ETPtype, viewsNum1, viewsNum2, X, Y, Z, RX, RY, RZ, Rdiff,
#                           density, RepErr_v1, RepErr_v2_e1, RepErr_v2_e2, RepErr_v3_e1, RepErr_v3_e2,
#                           residual_v1, residual_v2,TriAng_Max, TriAng_Avg, TriAng_Min, EpoMatch_Num,
#                           RE, RU, PA, IC, stdP, stdE, stdM, stdF, EF]
#     '''
#     crs = chunk.crs
#     M, T, R = MsTrans.getCoordTransMat(chunk)
#     ICTPsData = []
#     for i, point_id in enumerate(ICTPs_IdList):
#         Point = Points[point_id]
#         EpochTrack = ICTPsTracks[i]
#         Point3D_Retri = CTPsRetriCoord[i]
#         track_id = Point[-1]
#         ETPtype = CTPsTyp[i]
#         viewsNum1 = len(EpochTrack[0])
#         viewsNum2 = len(EpochTrack[1])
#         [X, Y, Z] = Point[0][1]  # in PJCS
#         [RX, RY, RZ] = MsTrans.transPointCoord_CLCS2PJCS(Point3D_Retri, crs, M)
#         Rdiff = np.linalg.norm(np.asarray([X - RX, Y - RY, Z - RZ]))
#         density = CTPsDen[i][1]
#         [point_id, RepErr_v1, RepErr_v2_e1, RepErr_v2_e2, RepErr_v3_e1, RepErr_v3_e2,
#          residual_v1, residual_v2_X, residual_v2_Z, residual_v2_Y, residual_v2_P, residual_v2_T,
#          TriAng_Max, TriAng_Avg, TriAng_Min, EpoMatch_Num] = CTPsQua[i][1:]
#         residual_v2_X = residual_v2_X * M.scale()
#         residual_v2_Z = residual_v2_Z * M.scale()
#         residual_v2_Y = residual_v2_Y * M.scale()
#         residual_v2_P = residual_v2_P * M.scale()
#         residual_v2_T = residual_v2_T * M.scale()
#         MCTPData = [point_id, track_id, ETPtype, viewsNum1, viewsNum2, X, Y, Z, RX, RY, RZ, Rdiff,
#                     density, RepErr_v1, RepErr_v2_e1, RepErr_v2_e2, RepErr_v3_e1, RepErr_v3_e2,
#                     residual_v1, residual_v2_X, residual_v2_Z, residual_v2_Y, residual_v2_P, residual_v2_T,
#                     TriAng_Max, TriAng_Avg, TriAng_Min, EpoMatch_Num]
#         ICTPsData.append(MCTPData)
#     return ICTPsData


def analyseICTPsNum(Points, point_ids, TracksEpoch):
    '''
    input:
        Points = {point_id:Point}
            ,where Point = [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, Track_id]
            ,where Coordinates=[[X, Y, Z],[X, Y, Z]]   # in CLCS & PJCS respectively
            ,where Covariance = 3x3 covariance matrix
            ,where StandardDeviation = [stdX, stdY, stdZ]  # in mm
            ,where Quality = [RE, RU, PA, IC]
        TracksEpoch = [[e1,e2,e3,e4], [e1,e2,e3,e4],…, [e1,e2,e3,e4]]
    output:
        CTPs_Number = [numCommonTiePoints, numRetained_Both, numRetained_1, numRetained_2, numRemoved, numEpoch1, numEpoch2]
    Note that this function only used for pairwise co-align.
    '''
    numRetained_Both = 0
    numRetained_1 = 0
    numRetained_2 = 0
    numRemoved = 0
    numEpoch1 = 0
    numEpoch2 = 0
    for track_id, epochNum in enumerate(TracksEpoch):
        # check
        if len(epochNum) != 2:
            raise Exception("[Script]    this script could only deal with pairwise co-align!")
        # skip invalid Track
        point_id = point_ids[track_id]
        if point_id == -1 or Points[point_id][5] == False:
            continue
        # judge which type dose this track belonging based on views number of each epoch
        e1 = epochNum[0]
        e2 = epochNum[1]
        if e1 > 1 and e2 > 1:
            numRetained_Both += 1
        elif e1 > 1 and e2 == 1:
            numRetained_1 += 1
        elif e1 == 1 and e2 > 1:
            numRetained_2 += 1
        elif e1 == 1 and e2 == 1:
            numRemoved += 1
        elif e1 == 0:
            numEpoch2 += 1
        elif e2 == 0:
            numEpoch1 += 1
        else:
            raise Exception("[Script]    the Track contains no points!")
    # total number of Common Tie Points
    numCommonTiePoints = numRetained_Both + numRetained_1 + numRetained_2 + numRemoved
    # numbers of each type points
    CTPs_Number = [numCommonTiePoints, numRetained_Both, numRetained_1, numRetained_2, numRemoved, numEpoch1, numEpoch2]
    return CTPs_Number


def analyseICTPsSpa(Points, weights, ICTPs_IdList):
    '''
    input:
        Points = {point_id:Point}
            ,where Point = [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, Track_id]
            ,where Coordinates=[[X, Y, Z],[X, Y, Z]]   # in CLCS & PJCS respectively
            ,where Covariance = 3x3 covariance matrix
            ,where StandardDeviation = [stdX, stdY, stdZ]  # in mm
            ,where Quality = [RE, RU, PA, IC]
        weights = [0, 0, 0, 0, 0, 0, 2, 4, 8, 16]
        ICTPs_IdList = [point_id, point_id,…,point_id]
    output:
        CTPs_Spatial = [score, numLevels, scoreLevels]
            ,where score = the total score of spacial distribution of CTPs
            ,where numLevels = [cells, cells, ...,cells],the number of cells of each level
            ,where scoreLevels = [score, score, ..., score],the score of each level
    This script provide functions for assessing the uniformity of spacial distribution by calcuate
    the spacial distribution score,SDS.
    Note that:
        the function need a third library: numpy, so the script use these function should be run in CMD
        X direct to East, Y direct to North
        the boundary and points used here are both in PJCS
    '''
    # step1  creat Pyramid
    Boundary_CLCS, Boundary_PJCS = Func_SpatialDistribution.getBoundary(Points)
    boundary_X_MAX = Boundary_PJCS[0]
    boundary_X_MIN = Boundary_PJCS[1]
    boundary_Y_MAX = Boundary_PJCS[2]
    boundary_Y_MIN = Boundary_PJCS[3]
    length_X = boundary_X_MAX - boundary_X_MIN
    length_Y = boundary_Y_MAX - boundary_Y_MIN

    # step2  calculate the matrixs
    matVariables = []
    numLevels = []
    scoreLevels = []
    for i in range(len(weights)):
        rows = math.ceil(length_Y / (2 ** i))
        cols = math.ceil(length_X / (2 ** i))
        if rows < 3 or cols < 3:
            break
        print('[Script]    level_{}:'.format(i), rows, cols)
        locals()['level_' + str(i)] = np.mat(np.zeros([rows, cols], dtype=np.uintc, order='C'))
        for point_id in ICTPs_IdList:
            Coordinates_PJCS = Points[point_id][0][1]
            X = Coordinates_PJCS[0]
            Y = Coordinates_PJCS[1]
            r = int((boundary_Y_MAX - float(Y)) / (2 ** i))  # row distance to origin(top-left)
            c = int((float(X) - boundary_X_MIN) / (2 ** i))  # col distance to origin(top-left)
            locals()['level_' + str(i)][r, c] = 1
        matVariables.append(locals()['level_' + str(i)])
    print('[Script]    Matrix Variables number: {}'.format(len(matVariables)))

    # step3  calculate SDS
    score = 0
    for i in range(len(matVariables)):
        cells = np.sum(matVariables[i])
        score += cells * weights[i]
        numLevels.append(cells)
        scoreLevels.append(cells * weights[i])
        print('[Script]    level_{}-cells: {}'.format(i, cells))
    CTPs_Spatial = [score, numLevels, scoreLevels]
    return CTPs_Spatial


def analyseICTPsType(ICTPs_IdList, Points, TracksEpoch):
    '''
    input:
        ICTPs_IdList = [point_id,point_id,...,point_id]
        Points
        TracksEpoch = [[e1,e2,e3,e4], [e1,e2,e3,e4],…, [e1,e2,e3,e4]] by track_id
    output:
        CTPs_Type = [ETP_Type,ETP_Type,...,ETP_Type]
    Note that this function only used for pairwise co-align.
    '''
    CTPs_Type = []
    for i, point_id in enumerate(ICTPs_IdList):
        TrackEpoch = TracksEpoch[Points[point_id][-1]]
        # this function only used for pairwise co-alignment!"
        # judge which type dose this track belonging based on views number of each epoch
        e1 = TrackEpoch[0]
        e2 = TrackEpoch[1]
        if e1 > 1 and e2 > 1:
            CTPs_Type.append('BothRetained')
        elif e1 > 1 and e2 == 1:
            CTPs_Type.append('Retained1')
        elif e1 == 1 and e2 > 1:
            CTPs_Type.append('Retained2')
        elif e1 == 1 and e2 == 1:
            CTPs_Type.append('Removed')
        else:
            raise Exception("[Script]    this epoch contains no views!")
    return CTPs_Type


def analyseICTPsDen(Points, ICTPs_IdList):
    '''
    input:
        Points = {point_id:Point}
            ,where Point = [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, Track_id]
            ,where Coordinates=[[X, Y, Z],[X, Y, Z]]   # in CLCS & PJCS respectively
            ,where Covariance = 3x3 covariance matrix
            ,where StandardDeviation = [stdX, stdY, stdZ]  # in mm
            ,where Quality = [RE, RU, PA, IC]
        ICTPs_IdList = [point_id, point_id,…,point_id]
    output:
        CTPs_Density = [ETPDensity, ETPDensity, ..., ETPDensity]
            ,where ETPDensity = [point_id, density]
     CTPsCoord should use the CLCS.
    '''
    ICTPsCoord_CLCS, ICTPsCoord_PJCS = Func_CommonTiePoints.getICTPsCoord(Points, ICTPs_IdList)
    CTPsCoord = np.asarray(ICTPsCoord_CLCS)
    values = np.vstack([CTPsCoord[:, 0], CTPsCoord[:, 1]])  # CTPsCoordX,CTPsCoordY
    densitys = gaussian_kde(values)(values)
    CTPs_Density = []
    for i, point_id in enumerate(ICTPs_IdList):
        density = densitys[i]
        CTPs_Density.append([point_id, density])
    return CTPs_Density


def analyseICTPsQua(chunk, ICTPsTracks, ICTPs_IdList, Cameras, Sensors):
    '''
    input:
        ICTPsTracks = [EpochTrack, EpochTrack,…, EpochTrack] (with index accompanied by ICTPsTracks_ids)
            ,where EpochTrack = [views_e1,views_e2,views_e3,views_e4],
            ,where views_e1 = [view,view,...,view]
            ,where view = [camera_id,projection_info],
            ,where projection_info = [projection.coord[0], projection.coord[1], projection.size, projection.track_id].
        ICTPs_IdList = [point_id,point_id,...,point_id]

        Tracks = [Track,Track,...,Track], a list that return Track by given track_id
            ,where Track = [view,view,...view],
            ,where view = [camera_id,projection_info]
            ,where projection_info = [projection.coord[0], projection.coord[1], projection.size, projection.track_id]
        Cameras
        Sensors
        TracksEpoch = [[e1,e2,e3,e4], [e1,e2,e3,e4],…, [e1,e2,e3,e4]]
        ICTPs_IdList = [point_id, point_id,…,point_id]
    output:
        CTPsQua = [CTPQua,CTPQua,...,CTPQua]
            ,where CTPQua = [point_id,
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
        CTPsRetriCoord = [Point3D,Point3D,...,Point3D]
            ,where Point3D = ndarray([X,Y,Z,1])
    Retriangulate CTPs and check its quality indexes
    '''
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)
    cameras = chunk.cameras

    CTPsQuality = []
    CTPsRetriCoord = []
    for EpochTrack_id, EpochTrack in enumerate(ICTPsTracks):
        # ['point_id', 'reprojection_error', 'triangulation_angle', 'epoch_match_number']
        point_id = ICTPs_IdList[EpochTrack_id]

        # reprojection_error version1
        views_all = [view for views in EpochTrack for view in views]

        Point3D_CLCS = Tri.Triangulate(views_all, Cameras, Sensors)
        error_pixel = Tri.PointRepErrSigma(Point3D_CLCS, views_all, Cameras, Sensors)

        # reprojection_error version2 and reprojection_error version3
        # note: this script only used for pairwise coalignment
        views_e1 = EpochTrack[0]
        views_e2 = EpochTrack[1]
        e1 = len(views_e1)
        e2 = len(views_e2)

        if e1 >= 2 and e2 >= 2:
            Point3D_e1 = Tri.Triangulate(views_e1, Cameras, Sensors)
            Point3D_e2 = Tri.Triangulate(views_e2, Cameras, Sensors)
            Point3D_e1_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e1, crs, M)
            Point3D_e2_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e2, crs, M)
            RepErr_v2_e1 = Tri.PointRepErrSigma(Point3D_e1, views_e1, Cameras, Sensors)
            RepErr_v2_e2 = Tri.PointRepErrSigma(Point3D_e2, views_e2, Cameras, Sensors)
            RepErr_v3_e1 = Tri.PointRepErrSigma(Point3D_e1, views_e2, Cameras, Sensors)
            RepErr_v3_e2 = Tri.PointRepErrSigma(Point3D_e2, views_e1, Cameras, Sensors)
            residual_X = Point3D_e1_PJCS[0] - Point3D_e2_PJCS[0]
            residual_Y = Point3D_e1_PJCS[1] - Point3D_e2_PJCS[1]
            residual_Z = Point3D_e1_PJCS[2] - Point3D_e2_PJCS[2]
            residual_P = np.linalg.norm(Point3D_e1_PJCS[:-2] - Point3D_e2_PJCS[:-2])
            residual_T = np.linalg.norm(Point3D_e1_PJCS[:-1] - Point3D_e2_PJCS[:-1])
        else:
            residual_X = 0
            residual_Y = 0
            residual_Z = 0
            residual_P = 0
            residual_T = 0
            if e1 >= 2 and e2 == 1:
                Point3D_e1 = Tri.Triangulate(views_e1, Cameras, Sensors)
                Point3D_e1_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e1, crs, M)
                Point3D_e2_PJCS = np.asarray([0, 0, 0, 0])
                RepErr_v2_e1 = Tri.PointRepErrSigma(Point3D_e1, views_e1, Cameras, Sensors)
                RepErr_v3_e1 = Tri.PointRepErrSigma(Point3D_e1, views_e2, Cameras, Sensors)
                RepErr_v2_e2 = 0
                RepErr_v3_e2 = 0
            elif e1 == 1 and e2 >= 2:
                Point3D_e2 = Tri.Triangulate(views_e2, Cameras, Sensors)
                Point3D_e1_PJCS = np.asarray([0, 0, 0, 0])
                Point3D_e2_PJCS = MsTrans.transPointCoord_CLCS2PJCS(Point3D_e2, crs, M)
                RepErr_v2_e2 = Tri.PointRepErrSigma(Point3D_e2, views_e2, Cameras, Sensors)
                RepErr_v3_e2 = Tri.PointRepErrSigma(Point3D_e2, views_e1, Cameras, Sensors)
                RepErr_v2_e1 = 0
                RepErr_v3_e1 = 0
            else:
                Point3D_e1_PJCS = np.asarray([0, 0, 0, 0])
                Point3D_e2_PJCS = np.asarray([0, 0, 0, 0])
                RepErr_v2_e1 = 0
                RepErr_v2_e2 = 0
                RepErr_v3_e1 = 0
                RepErr_v3_e2 = 0

        # triangulation_angle
        TriAng_List = Tri.PointTriAngList(Point3D_CLCS, views_all, Cameras)
        TriAng_Max = max(TriAng_List)
        TriAng_Avg = np.mean(TriAng_List)
        TriAng_Min = min(TriAng_List)

        # epoch_match_number
        EpoMatch_Num = len(views_e1) * len(views_e2)

        CTPQuality = [point_id,
                      Point3D_e1_PJCS[0], Point3D_e1_PJCS[1], Point3D_e1_PJCS[2],
                      Point3D_e2_PJCS[0], Point3D_e2_PJCS[1], Point3D_e2_PJCS[2],
                      residual_X, residual_Y, residual_Z, residual_P, residual_T,
                      error_pixel, RepErr_v3_e1, RepErr_v3_e2, RepErr_v2_e1, RepErr_v2_e2,
                      TriAng_Max, TriAng_Avg, TriAng_Min,
                      len(views_all), e1, e2, EpoMatch_Num]
        CTPsQuality.append(CTPQuality)
        CTPsRetriCoord.append(Point3D_CLCS)
    return CTPsQuality, CTPsRetriCoord


# def analyseICTPsQua2(chunk, ICTPs_IdList):
#     '''
#     input:
#         ICTPs_IdList = [point_id,point_id,...,point_id]
#     output:
#         CTPs_MetashapeQuality = [ETPMsQuality,ETPMsQuality,...,ETPMsQuality]
#             ,where ETPMsQuality = [point_id, RE, RU, PA, IC, stdP, stdE, stdM, stdF, EF]
#     Read CTPs quality indexes from Metashape
#     '''
#     CTPsMetashapeQuality = []
#     CTPsCriterion = MsPoint.getPointsCriterion(chunk, IdList=ICTPs_IdList)
#     CTPsCov = MsPoint.getPointsCov_analysis(chunk, IdList=ICTPs_IdList)
#     for i, point_id in enumerate(ICTPs_IdList):
#         [RE, RU, PA, IC] = CTPsCriterion[i]
#         [stdP, stdE, stdM, stdF, EF] = CTPsCov[i]
#         ETPMsQua = [point_id, RE, RU, PA, IC, stdP, stdE, stdM, stdF, EF]
#         CTPsMetashapeQuality.append(ETPMsQua)
#     return CTPsMetashapeQuality
# def exportICTPsQua(CTPs, path):
#     '''
#     input:
#         CTPs = [CTP,CTP,...,CTP]
#             ,where CTP = [point_id, track_id, ETPtype, viewsNum1, viewsNum2, X, Y, Z, RX, RY, RZ, Rdiff, density,
#                             0           1        2         3          4      5  6  7  8   9   10   11     12
#                           RepErr_v1, RepErr_v2_e1, RepErr_v2_e2, RepErr_v3_e1, RepErr_v3_e2,
#                             13             14           15            16            17
#                           residual_v1, residual_v2_X, residual_v2_Z, residual_v2_Y, residual_v2_P, residual_v2_T,
#                               18            19              20             21              22           23
#                           TriAng_Max, TriAng_Avg, TriAng_Min, EpoMatch_Num, RE, RU, PA, IC, stdP, stdE, stdM, stdF, EF]
#                               24          25           26           27      28  29  30  31   32    33    34    35   36
#     output:
#         _CTPsQuality.txt
#     '''
#     reportFile = open(path, "w")
#     fwriter = csv.writer(reportFile, delimiter='\t', lineterminator='\n')
#     fwriter.writerow(
#         ['point_id', 'track_id', 'ETPtype', 'viewsNum1', 'viewsNum2', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'Rdiff',
#          'density', 'RepErr_v1', 'RepErr_v2_e1', 'RepErr_v2_e2', 'RepErr_v3_e1', 'RepErr_v3_e2',
#          'residual_v1', 'residual_v2_X, residual_v2_Z, residual_v2_Y, residual_v2_P, residual_v2_T',
#          'TriAng_Max', 'TriAng_Avg', 'TriAng_Min', 'EpoMatch_Num', 'RE', 'RU', 'PA', 'IC', 'stdP',
#          'stdE', 'stdM', 'stdF', 'EF'])
#     for CTP in CTPs:
#         fwriter.writerow(['{}'.format(CTP[0]), '{}'.format(CTP[1]), '{}'.format(CTP[2]),
#                           '{}'.format(CTP[3]), '{}'.format(CTP[4]),
#                           '{0:0.5f}'.format(CTP[5]), '{0:0.5f}'.format(CTP[6]), '{0:0.5f}'.format(CTP[7]),
#                           '{0:0.5f}'.format(CTP[8]), '{0:0.5f}'.format(CTP[9]), '{0:0.5f}'.format(CTP[10]),
#                           '{0:0.5f}'.format(CTP[11]), '{0:0.5f}'.format(CTP[12]),
#                           '{0:0.5f}'.format(CTP[13]), '{0:0.5f}'.format(CTP[14]), '{0:0.5f}'.format(CTP[15]),
#                           '{0:0.5f}'.format(CTP[16]), '{0:0.5f}'.format(CTP[17]),
#                           '{0:0.5f}'.format(CTP[18]), '{0:0.5f}'.format(CTP[19]), '{0:0.5f}'.format(CTP[20]),
#                           '{0:0.5f}'.format(CTP[21]), '{0:0.5f}'.format(CTP[22]), '{0:0.5f}'.format(CTP[23]),
#                           '{0:0.3f}'.format(CTP[24]), '{0:0.3f}'.format(CTP[25]), '{0:0.3f}'.format(CTP[26]),
#                           '{}'.format(CTP[27]), '{0:0.5f}'.format(CTP[28]), '{0:0.5f}'.format(CTP[29]),
#                           '{0:0.5f}'.format(CTP[30]), '{}'.format(CTP[31]), '{0:0.5f}'.format(CTP[32]),
#                           '{0:0.5f}'.format(CTP[32]), '{0:0.5f}'.format(CTP[33]), '{0:0.5f}'.format(CTP[34]),
#                           '{0:0.5f}'.format(CTP[35])])
#     reportFile.close()

def reportICTPsQua(CTPsQua, path):
    '''
    input:
        CTPsQua = [CTPQua,CTPQua,...,CTPQua]
            ,where CTPQua = [point_id,
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
        _CTPsQua.txt
    '''
    # get statistic of residuals
    residual_X_List = [AllQua[7] for AllQua in CTPsQua if AllQua[11] != 0]
    residual_Y_List = [AllQua[8] for AllQua in CTPsQua if AllQua[11] != 0]
    residual_Z_List = [AllQua[9] for AllQua in CTPsQua if AllQua[11] != 0]
    residual_P_List = [AllQua[10] for AllQua in CTPsQua if AllQua[11] != 0]
    residual_T_List = [AllQua[11] for AllQua in CTPsQua if AllQua[11] != 0]
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
    Distrib_ErrorPixel = FcSta.listStatistic([AllQua[12] for AllQua in CTPsQua])
    Distrib_RepErr_v3_e1 = FcSta.listStatistic([AllQua[13] for AllQua in CTPsQua if AllQua[11] != 0])
    Distrib_RepErr_v3_e2 = FcSta.listStatistic([AllQua[14] for AllQua in CTPsQua if AllQua[11] != 0])
    Distrib_RepErr_v2_e1 = FcSta.listStatistic([AllQua[15] for AllQua in CTPsQua if AllQua[15] != 0])
    Distrib_RepErr_v2_e2 = FcSta.listStatistic([AllQua[16] for AllQua in CTPsQua if AllQua[16] != 0])
    Distrib_TriAng_Max = FcSta.listStatistic([AllQua[17] for AllQua in CTPsQua])
    Distrib_TriAng_Avg = FcSta.listStatistic([AllQua[18] for AllQua in CTPsQua])
    Distrib_TriAng_Min = FcSta.listStatistic([AllQua[19] for AllQua in CTPsQua])
    Distrib_ViewNum = FcSta.listStatistic([AllQua[20] for AllQua in CTPsQua])
    Distrib_ViewNum1 = FcSta.listStatistic([AllQua[21] for AllQua in CTPsQua])
    Distrib_ViewNum2 = FcSta.listStatistic([AllQua[22] for AllQua in CTPsQua])
    Distrib_EpoMatch_Num = FcSta.listStatistic([AllQua[23] for AllQua in CTPsQua])

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


def reportICTPsInfo(CTPsNum, CTPsSpa, weights, path):
    '''
    input:
        CTPsNum = [numCommonTiePoints, numRetained_Both, numRetained_1, numRetained_2, numRemoved, numEpoch1, numEpoch2]
        CTPsSpa = [score, numLevels, scoreLevels]
            ,where score = the total score of spacial distribution of CTPs
            ,where numLevels = [cells, cells, ...,cells],the number of cells of each level
            ,where scoreLevels = [score, score, ..., score],the score of each level
        weights = [0, 0, 0, 0, 0, 0, 2, 4, 8, 16]
    output:
        _ICTPsInfo.txt
    '''
    # write file
    reportFile = open(path, "w")
    fwriter = csv.writer(reportFile, delimiter='\t', lineterminator='\n')
    # number of CTPs
    fwriter.writerow(['numPoints', CTPsNum[0] + CTPsNum[-2] + CTPsNum[-1]])
    fwriter.writerow(['numCommonTiePoints', CTPsNum[0]])
    fwriter.writerow(['numRetained_Both', CTPsNum[1]])
    fwriter.writerow(['numRetained_1', CTPsNum[2]])
    fwriter.writerow(['numRetained_2', CTPsNum[3]])
    fwriter.writerow(['numRemoved', CTPsNum[4]])
    fwriter.writerow(['numEpoch1', CTPsNum[5]])
    fwriter.writerow(['numEpoch2', CTPsNum[6]])
    # spacial distribution of CTPs
    fwriter.writerow(['spacial distribution score', CTPsSpa[0]])
    fwriter.writerow(weights)
    fwriter.writerow(CTPsSpa[1])
    fwriter.writerow(CTPsSpa[2])
    reportFile.close()
    return None


def deleteICTPs(ICTPs_IdList, chunk, Points, TracksEpoch, ICTPsTypes):
    '''
    input:
        chunk
        ICTPs_IdList
        TracksEpoch = [[e1,e2,e3,e4], [e1,e2,e3,e4],…, [e1,e2,e3,e4]]
        ICTPsTypes = ['BothRetained', 'Retained1','Retained2','Removed']
    output:
        delete the selected points in chunk.
    '''
    CTPs_Type = analyseICTPsType(ICTPs_IdList, Points, TracksEpoch)
    if ICTPsTypes:
        for i, point_id in enumerate(ICTPs_IdList):
            if CTPs_Type[i] in ICTPsTypes:
                chunk.point_cloud.points[point_id].valid = False
    else:
        for i, point_id in enumerate(ICTPs_IdList):
            chunk.point_cloud.points[point_id].valid = False

