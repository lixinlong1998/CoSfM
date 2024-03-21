from config import *
import csv
import numpy as np
from src_CFTM import ConnectData_Metashape
from src_CFTM import Func_CommonTiePoints
from src_CFTM import ConnectData_COLMAP as Colmap
from src_CFTM import Func_Files as FcFile
from src_Metashape import FuncMs_Analysis
from src_Metashape import FuncMs_Marker as MsMarker

'''
path_Matches:
    [MatchesReport_cube_e1, MatchesReport_cube_e2]
    MatchesReport_cube = [MatchesReport_cube_pair_id,...]
    MatchesReport_cube_pair_id = [(camera_id1, camera_id2), Matches, inlierMatches]
        ,where Matches = np.array([[keypoint_id, keypoint_id],…])
        ,where inlierMatches = np.array([[keypoint_id, keypoint_id],…])
path_TracksReport:
    [TracksReport_cube_e1, TracksReport_cube_e2]
    TracksReport_cube = [TrackReport,...]
        ,where TrackReport = [views,edges]
        ,where views = np.array([[camera_id, keypoint_id],…])
        ,where edges = np.array([[view_id, view_id],…])
path_FilterTracksReport:
    [FilterTracksReport_cube_e1, FilterTracksReport_cube_e2]
    TracksReport_cube = [TrackReport,...]
        ,where TrackReport = [views,edges]
        ,where views = np.array([[camera_id, keypoint_id],…])
        ,where edges = np.array([[view_id, view_id],…])
path_CommonTracksReport:
    CommonTracksMatches = [CommonTrackMatches,...]
        ,where CommonTrackMatches = [CrossMatch,...]
        ,where CrossMatch = [feature1, feature2, MatchQuality]
        ,where feature = [camera_id, keypoint_id]
        ,where MatchQuality = [distance,angles],distance is the L2 normal distance between descriptors of match.
path_CTPs:
    for i, CommonTrack in enumerate(CommonTracks):
        fwriter.writerow([CommonTrack, CommonTracksMatches[i]])
    CommonTrack = [view,view,...,view]
        ,where view = [camera_id, projection_info, index_info]
        ,where projection_info = [u, v]
        ,where index_info = [identifier, keypoints_id, originTrack_id]
    CommonTrackMatches = [CrossMatch,…]
        ,where CrossMatch = [feature1, feature2, MatchQuality]
        ,where feature = [camera_id, keypoint_id]
        ,where MatchQuality = [distance,angles],distance is the L2 normal distance between descriptors of match.
'''


def analysisCTPs(chunk, args):
    starttime0 = time.perf_counter()
    print('[Script]    Analysis CTPs...')

    # [0] unpack arguments
    project_name = args["project_name"]
    data_package_path = args["data_package_path"]
    report_path = args["report_path"]
    os.makedirs(report_path, exist_ok=True)
    MarkerGroupName = 'CTPs'

    # [1]  prepare data
    camera_ids = FcFile.readDictionary_Int_Int(data_package_path + '/camera_ids.txt')
    CamerasEpoch = FcFile.readCamerasEpoch(data_package_path + '/CamerasEpoch.txt')

    # [2]  read GCTPs database
    GCTPsData = MsMarker.getMarkersData_Analyse(chunk, MarkerGroupName, mode='uncheck')

    # [3]  Analyse GCTPs
    GCTPsQua = MsMarker.getMarkersQua_Analyse(chunk, GCTPsData, camera_ids, CamerasEpoch)
    # CTPsAllQua = MsMarker.getMarkersAllQua(chunk, CTPsData, camera_ids, CamerasEpoch, 0,
    #                                        MarkerFormat='Analyse', Triangulation='linear')

    # [4]  generate report and export results
    path_GCTPsQua = os.path.join(report_path, os.path.splitext(project_name)[0] + '_GCTPsQua.txt')
    path_GCTPsQuaSta = os.path.join(report_path, os.path.splitext(project_name)[0] + '_GCTPsQuaSta.txt')
    MsMarker.exportMarkersQua(GCTPsQua, path_GCTPsQua)
    MsMarker.reportMarkersQua(GCTPsQua, path_GCTPsQuaSta)
    print('[Script][TimeCost]    analysis  analysed:', time.perf_counter() - starttime0)


def analysisCPs(chunk, args):
    starttime0 = time.perf_counter()
    print('[Script]    Analyse check points...')

    # [0] unpack arguments
    project_name = args["project_name"]
    check_points_path = args["check_points_path"]
    data_package_path = args["data_package_path"]
    report_path = args["report_path"]
    os.makedirs(report_path, exist_ok=True)

    # [1]  prepare data
    camera_ids = FcFile.readDictionary_Int_Int(data_package_path + '/camera_ids.txt')
    CamerasEpoch = FcFile.readCamerasEpoch(data_package_path + '/CamerasEpoch.txt')

    # [2]  Analyse check points
    CPsData = MsMarker.importMarkersData_Analyse(check_points_path, MarkerList=[])
    CPsQua = MsMarker.getMarkersQua_Analyse(chunk, CPsData, camera_ids, CamerasEpoch)

    # [3]  Generate report and export results
    path_CPsQua = os.path.join(report_path, os.path.splitext(project_name)[0] + '_CPsQua.txt')
    path_CPsQuaSta = os.path.join(report_path, os.path.splitext(project_name)[0] + '_CPsQuaSta.txt')
    MsMarker.exportMarkersQua(CPsQua, path_CPsQua)
    MsMarker.reportMarkersQua(CPsQua, path_CPsQuaSta)


def analysisMatches(args):
    starttime0 = time.perf_counter()
    print('[Script]    Analysis Matches...')

    # [0] unpack arguments
    CTP_path = args["CTP_path"]
    report_path = args["report_path"]
    os.makedirs(report_path, exist_ok=True)

    # [1] Import data
    MatchesReport, TracksReport, FilterTracksReport, CommonTracksReport = FuncMs_Analysis.loadFolder_CTPs(CTP_path)

    # [2] Parse data
    # Read inner matches: store as pairs containing original and inlier matches
    MatchesReport_e1, MatchesReport_e2 = FuncMs_Analysis.readInnerMatches(MatchesReport)
    # Inner matches: store as original and inlier matches: original matches store e1e2, inlier matches store e1e2
    MatchesSta, InlierMatchesSta = FuncMs_Analysis.reshapeMatchesStructure(MatchesReport_e1, MatchesReport_e2)
    # Read cross-temporal matches: convert from track-based storage to pair-based storage
    CrossMatchesReport_dict = FuncMs_Analysis.readCrossMatches(CommonTracksReport)
    # Read the number of cross-temporal matches for common tracks
    CTracksCrMatchReport = FuncMs_Analysis.readCommonTracksCrossMatch(CommonTracksReport)
    # Read inner tracks
    InnerTracksReport_e1, InnerTracksReport_e2 = FuncMs_Analysis.readInnerTracks(TracksReport)
    InnerFiltTracksReport_e1, InnerFiltTracksReport_e2 = FuncMs_Analysis.readInnerTracks(FilterTracksReport)

    # [3] Analyze data
    '''
    # //////Inner pair count
    # //////Total inner matches count (original/inlier)
    # //////Average inner matches count per pair (original/inlier)
    # //////Average inlier ratio per inner pair
    # //////Cross-temporal pair count
    # //////Total cross-temporal matches count (count of cross-temporal matches with CFTM enabled)
    # //////Percentage of cross-temporal matches in total (enabled) matches count
    # //////Average cross-temporal matches count per pair
    '''
    # Inner pair count (original/inlier)
    Epoch1_PairNum = [len(MatchesSta[0]), len(InlierMatchesSta[0])]
    Epoch2_PairNum = [len(MatchesSta[1]), len(InlierMatchesSta[1])]
    # Total inner matches count (original/inlier)
    Epoch1_MatchNum = [sum(MatchesSta[0]), sum(InlierMatchesSta[0])]
    Epoch2_MatchNum = [sum(MatchesSta[1]), sum(InlierMatchesSta[1])]
    # Average inner matches count per pair (original/inlier)
    Epoch1_PairMatcNum_AVG = Epoch1_MatchNum[0] / Epoch1_PairNum[0]
    Epoch2_PairMatcNum_AVG = Epoch2_MatchNum[0] / Epoch2_PairNum[0]
    # Average inlier ratio per inner pair
    Epoch1_PairInlierMatcRatio_AVG = Epoch1_MatchNum[1] / Epoch1_MatchNum[0]
    Epoch2_PairInlierMatcRatio_AVG = Epoch2_MatchNum[1] / Epoch2_MatchNum[0]
    # Cross-temporal pair count
    Common_PairNum = len(CrossMatchesReport_dict)
    # Total cross-temporal matches count
    Common_MatchNum1 = sum([len(CrossMatches) for pair_id, CrossMatches in CrossMatchesReport_dict.items()])
    Common_MatchNum2 = sum(CTracksCrMatchReport)
    assert Common_MatchNum1 == Common_MatchNum2
    # Percentage of cross-temporal matches in total (enabled) matches count
    Common_MatchPortion = Common_MatchNum1 / (Common_MatchNum1 + Epoch1_MatchNum[1] + Epoch2_MatchNum[1])
    # Average cross-temporal matches count per pair
    Common_PairMatcNum_AVG = Common_MatchNum1 / Common_PairNum
    '''
    # //////Inner track count (before RANSAC)
    # //////Inner track count (after RANSAC)
    # //////Cross-temporal track count
    # //////Average cross-temporal matches count per cross-temporal track
    # //////Histogram of cross-temporal matches count for cross-temporal tracks: hist(CTracksCrMatchReport)
    '''
    # Inner track count (before RANSAC)
    Epoch1_TrackNum = len(InnerTracksReport_e1)
    Epoch2_TrackNum = len(InnerTracksReport_e2)
    # Inner track count (after RANSAC)
    Epoch1_InlierTrackNum = len(InnerFiltTracksReport_e1)
    Epoch2_InlierTrackNum = len(InnerFiltTracksReport_e2)
    # Cross-temporal track count
    Common_TrackNum = len(CTracksCrMatchReport)
    # Average cross-temporal matches count per cross-temporal track
    Common_TrackCrosMatcNum_AVG = Common_MatchNum1 / len(CTracksCrMatchReport)
    # Histogram of cross-temporal matches count for cross-temporal tracks: hist(CTracksCrMatchReport)

    # [4] Display results
    print('Epoch1_PairNum:', Epoch1_PairNum)
    print('Epoch2_PairNum:', Epoch2_PairNum)
    print('Epoch1_MatchNum:', Epoch1_MatchNum)
    print('Epoch2_MatchNum:', Epoch2_MatchNum)
    print('Epoch1_Epoch1_PairMatcNum_AVG:', Epoch1_PairMatcNum_AVG)
    print('Epoch2_Epoch2_PairMatcNum_AVG:', Epoch2_PairMatcNum_AVG)
    print('Epoch1_InlierMatchRatio:{0:0.3f}%'.format(Epoch1_PairInlierMatcRatio_AVG * 100))
    print('Epoch2_InlierMatchRatio:{0:0.3f}%'.format(Epoch2_PairInlierMatcRatio_AVG * 100))
    print('Common_PairNum:', Common_PairNum)
    print('Common_MatchNum1:', Common_MatchNum1)
    print('Common_MatchNum2:', Common_MatchNum2)
    print('Common_MatchPortion:', Common_MatchPortion)
    print('Common_PairMatcNum_AVG:', Common_PairMatcNum_AVG)
    print('Epoch1_TrackNum:', Epoch1_TrackNum)
    print('Epoch2_TrackNum:', Epoch2_TrackNum)
    print('Epoch1_InlierTrackNum:', Epoch1_InlierTrackNum)
    print('Epoch2_InlierTrackNum:', Epoch2_InlierTrackNum)
    print('Common_TrackNum:', Common_TrackNum)
    print('Common_TrackCrosMatcNum_AVG:', Common_TrackCrosMatcNum_AVG)

    # [5] Output results to npy
    Analyse_path = report_path + "/MatchesReport_CFTM"
    CFTMAnalyseReport = [[Epoch1_PairNum, Epoch2_PairNum, Epoch1_MatchNum, Epoch2_MatchNum],
                         [Epoch1_PairMatcNum_AVG, Epoch2_PairMatcNum_AVG],
                         [Epoch1_PairInlierMatcRatio_AVG, Epoch2_PairInlierMatcRatio_AVG],
                         [Common_PairNum, Common_MatchNum1, Common_MatchNum2, Common_MatchPortion,
                          Common_PairMatcNum_AVG],
                         [Epoch1_TrackNum, Epoch2_TrackNum, Epoch1_InlierTrackNum, Epoch2_InlierTrackNum],
                         [Common_TrackNum, Common_TrackCrosMatcNum_AVG]]
    np.save(Analyse_path + '.npy', CFTMAnalyseReport)

    # [6] Output results to txt
    File = open(Analyse_path + '.txt', "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    fwriter.writerow([Epoch1_PairNum, Epoch2_PairNum, Epoch1_MatchNum, Epoch2_MatchNum])
    fwriter.writerow([Epoch1_PairMatcNum_AVG, Epoch2_PairMatcNum_AVG])
    fwriter.writerow([Epoch1_PairInlierMatcRatio_AVG, Epoch2_PairInlierMatcRatio_AVG])
    fwriter.writerow([Common_PairNum, Common_MatchNum1, Common_MatchNum2, Common_MatchPortion, Common_PairMatcNum_AVG])
    fwriter.writerow([Epoch1_TrackNum, Epoch2_TrackNum, Epoch1_InlierTrackNum, Epoch2_InlierTrackNum])
    fwriter.writerow([Common_TrackNum, Common_TrackCrosMatcNum_AVG])
    File.close()

    print('[Script][TimeCost]    Matches analyzed:', time.perf_counter() - starttime0)


def compareMatches(chunk, args):
    starttime0 = time.perf_counter()
    print('[Script]    Compare Matches between co-alignment and CFTM...')

    # [0] unpack arguments
    data_package_path = args["data_package_path"]
    feature_path = args["feature_path"]
    report_path = args["report_path"]
    epoch_mode = args["epoch_mode"]
    os.makedirs(report_path, exist_ok=True)

    # [1]  prepare data
    camera_ids = FcFile.readDictionary_Int_Int(data_package_path + '/camera_ids.txt')
    cameraPaths = FcFile.readDictionary_Str_Int(data_package_path + '/cameraPaths.txt')
    points = chunk.point_cloud.points
    point_ids = ConnectData_Metashape.getPointIds(chunk)
    if epoch_mode == "DATE":
        epochs = Func_CommonTiePoints.getEpoch_byDATE(camera_ids)
    elif epoch_mode == "FOLDER":
        epochs = Func_CommonTiePoints.getEpoch_byFOLDER(cameraPaths)
    Tracks = ConnectData_Metashape.getTracks(chunk)
    Cameras = FcFile.readDictionary_Int_List(data_package_path + '/Cameras.txt')
    cameraPaths = FcFile.readDictionary_Str_Int(data_package_path + '/cameraPaths.txt')
    CamerasEpoch = FcFile.readCamerasEpoch(data_package_path + '/CamerasEpoch.txt')
    TracksEpoch = Func_CommonTiePoints.analyseTracks_Epoch(Tracks, epochs, CamerasEpoch)

    # [2] 在colmap中进行穷尽匹配以获得同名点信息；注意穷尽匹配会非常耗时
    chunkMatches = Colmap.extractMatch(Cameras, cameraPaths, feature_path)
    loaded_data_chunkDescriptors = np.load(os.path.join(data_package_path, 'chunkDescriptors.npz'))
    chunkDescriptors = loaded_data_chunkDescriptors['arr_0']
    loaded_data_chunkDescriptors.close()
    chunkInlierMatches, chunkGeometry = Colmap.extractInlierMatch(Cameras, cameraPaths, feature_path)

    # [3] 记录matches中每个匹配的距离
    MatchesSta, MatchesDis = FuncMs_Analysis.getMatchAnalyse(chunkMatches, chunkDescriptors, camera_ids, CamerasEpoch,
                                                             epochs)
    InlierMatchesSta, InlierMatchesDis = FuncMs_Analysis.getMatchAnalyse(chunkInlierMatches, chunkDescriptors,
                                                                         camera_ids, CamerasEpoch, epochs)
    MatchesReport_Epoch1, MatchesReport_Epoch2, MatchesReport_Common = FuncMs_Analysis.alignAllMatch(chunkMatches,
                                                                                                     chunkInlierMatches,
                                                                                                     camera_ids,
                                                                                                     CamerasEpoch,
                                                                                                     epochs)
    Epoch1_TrackReport, Epoch2_TrackReport, Common_TrackReport = FuncMs_Analysis.getTracksReport(chunk, TracksEpoch,
                                                                                                 points, point_ids,
                                                                                                 CriteriasThreshold=[
                                                                                                     0.5, 20, 8])

    # [4] 分析数据
    # 像对的数量
    Epoch1_PairNum = [len(MatchesSta[0]), len(InlierMatchesSta[0])]
    Epoch2_PairNum = [len(MatchesSta[1]), len(InlierMatchesSta[1])]
    Common_PairNum = [len(MatchesSta[2]), len(InlierMatchesSta[2])]
    # matches的数量
    Epoch1_MatchNum = [sum(MatchesSta[0]), sum(InlierMatchesSta[0])]
    Epoch2_MatchNum = [sum(MatchesSta[1]), sum(InlierMatchesSta[1])]
    Common_MatchNum = [sum(MatchesSta[2]), sum(InlierMatchesSta[2])]
    # 像对的平均原始匹配数
    Epoch1_PairMatchNum_AVG = Epoch1_MatchNum[0] / Epoch1_PairNum[0]
    Epoch2_PairMatchNum_AVG = Epoch2_MatchNum[0] / Epoch2_PairNum[0]
    Common_PairMatchNum_AVG = Common_MatchNum[0] / Common_PairNum[0]
    # 像对的平均过验匹配数(忽略内点匹配为0的像对)
    Epoch1_PairInlierMatchNum_AVG = Epoch1_MatchNum[1] / Epoch1_PairNum[1]
    Epoch2_PairInlierMatchNum_AVG = Epoch2_MatchNum[1] / Epoch2_PairNum[1]
    Common_PairInlierMatchNum_AVG = Common_MatchNum[1] / Common_PairNum[1]
    # 像对的加权平均匹配过验率
    Epoch1_PairInlierMatchRatio_WAVG = Epoch1_MatchNum[1] / Epoch1_MatchNum[0]
    Epoch2_PairInlierMatchRatio_WAVG = Epoch2_MatchNum[1] / Epoch2_MatchNum[0]
    Common_PairInlierMatchRatio_WAVG = Common_MatchNum[1] / Common_MatchNum[0]
    # 像对的平均匹配过验率(忽略内点匹配为0的像对)
    Epoch1_PairInlierMatchRatio_AVG = FuncMs_Analysis.getPairInlierMatchRatio(MatchesReport_Epoch1)
    Epoch2_PairInlierMatchRatio_AVG = FuncMs_Analysis.getPairInlierMatchRatio(MatchesReport_Epoch2)
    Common_PairInlierMatchRatio_AVG = FuncMs_Analysis.getPairInlierMatchRatio(MatchesReport_Common)
    # 各类原始匹配占比
    all_MatchNum = sum([sum(species) for species in MatchesSta])
    Epoch1_MatchPortion = Epoch1_MatchNum[0] / all_MatchNum
    Epoch2_MatchPortion = Epoch2_MatchNum[0] / all_MatchNum
    Common_MatchPortion = Common_MatchNum[0] / all_MatchNum
    # 各类过验匹配占比
    all_InlierMatchNum = sum([sum(species) for species in InlierMatchesSta])
    Epoch1_InlierMatchPortion = Epoch1_MatchNum[1] / all_InlierMatchNum
    Epoch2_InlierMatchPortion = Epoch2_MatchNum[1] / all_InlierMatchNum
    Common_InlierMatchPortion = Common_MatchNum[1] / all_InlierMatchNum
    # Track数量
    Epoch1_TrackNum = len(Epoch1_TrackReport)
    Epoch2_TrackNum = len(Epoch2_TrackReport)
    Common_TrackNum = len(Common_TrackReport)
    # valid track数量
    Epoch1_ValidTrackNum = sum([TR[3] for TR in Epoch1_TrackReport])
    Epoch2_ValidTrackNum = sum([TR[3] for TR in Epoch2_TrackReport])
    Common_ValidTrackNum = sum([TR[3] for TR in Common_TrackReport])
    # ER 后 valid track数量
    Epoch1_ERValidTrackNum = sum([TR[4] for TR in Epoch1_TrackReport])
    Epoch2_ERValidTrackNum = sum([TR[4] for TR in Epoch2_TrackReport])
    Common_ERValidTrackNum = sum([TR[4] for TR in Common_TrackReport])
    # valid track数量在所有track中的占比
    Epoch1_ValidTrackPortion = Epoch1_ValidTrackNum / Epoch1_TrackNum
    Epoch2_ValidTrackPortion = Epoch2_ValidTrackNum / Epoch2_TrackNum
    Common_ValidTrackPortion = Common_ValidTrackNum / Common_TrackNum
    # ER 后 valid track数量在所有track中的占比
    Epoch1_ERValidTrackPortion = Epoch1_ERValidTrackNum / Epoch1_TrackNum
    Epoch2_ERValidTrackPortion = Epoch2_ERValidTrackNum / Epoch2_TrackNum
    Common_ERValidTrackPortion = Common_ERValidTrackNum / Common_TrackNum
    print(Epoch1_ValidTrackPortion)
    print(Epoch2_ValidTrackPortion)
    print(Common_ValidTrackPortion)
    print(Epoch1_ERValidTrackPortion)
    print(Epoch2_ERValidTrackPortion)
    print(Common_ERValidTrackPortion)

    # [5] 输出结果
    Analyse_path = report_path + "/MatchesReport_Coalign"
    # to npy
    CoalignAnalyseReport = [[Epoch1_PairNum, Epoch2_PairNum, Common_PairNum],
                            [Epoch1_MatchNum, Epoch2_MatchNum, Common_MatchNum],
                            [Epoch1_PairMatchNum_AVG, Epoch2_PairMatchNum_AVG, Common_PairMatchNum_AVG],
                            [Epoch1_PairInlierMatchNum_AVG, Epoch2_PairInlierMatchNum_AVG,
                             Common_PairInlierMatchNum_AVG],
                            [Epoch1_PairInlierMatchRatio_WAVG, Epoch2_PairInlierMatchRatio_WAVG,
                             Common_PairInlierMatchRatio_WAVG],
                            [Epoch1_PairInlierMatchRatio_AVG, Epoch2_PairInlierMatchRatio_AVG,
                             Common_PairInlierMatchRatio_AVG],
                            [Epoch1_MatchPortion, Epoch2_MatchPortion, Common_MatchPortion],
                            [Epoch1_InlierMatchPortion, Epoch2_InlierMatchPortion, Common_InlierMatchPortion],
                            [Epoch1_TrackNum, Epoch2_TrackNum, Common_TrackNum],
                            [Epoch1_ValidTrackNum, Epoch2_ValidTrackNum, Common_ValidTrackNum],
                            [Epoch1_ERValidTrackNum, Epoch2_ERValidTrackNum, Common_ERValidTrackNum],
                            [Epoch1_ValidTrackPortion, Epoch2_ValidTrackPortion, Common_ValidTrackPortion],
                            [Epoch1_ERValidTrackPortion, Epoch2_ERValidTrackPortion, Common_ERValidTrackPortion]]
    np.save(Analyse_path + '.npy', CoalignAnalyseReport)

    # to txt
    File = open(Analyse_path + '.txt', "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    fwriter.writerow([Epoch1_PairNum, Epoch2_PairNum, Common_PairNum])
    fwriter.writerow([Epoch1_MatchNum, Epoch2_MatchNum, Common_MatchNum])
    fwriter.writerow([Epoch1_PairMatchNum_AVG, Epoch2_PairMatchNum_AVG, Common_PairMatchNum_AVG])
    fwriter.writerow([Epoch1_PairInlierMatchNum_AVG, Epoch2_PairInlierMatchNum_AVG, Common_PairInlierMatchNum_AVG])
    fwriter.writerow(
        [Epoch1_PairInlierMatchRatio_WAVG, Epoch2_PairInlierMatchRatio_WAVG, Common_PairInlierMatchRatio_WAVG])
    fwriter.writerow(
        [Epoch1_PairInlierMatchRatio_AVG, Epoch2_PairInlierMatchRatio_AVG, Common_PairInlierMatchRatio_AVG])
    fwriter.writerow([Epoch1_MatchPortion, Epoch2_MatchPortion, Common_MatchPortion])
    fwriter.writerow([Epoch1_InlierMatchPortion, Epoch2_InlierMatchPortion, Common_InlierMatchPortion])
    fwriter.writerow([Epoch1_TrackNum, Epoch2_TrackNum, Common_TrackNum])
    fwriter.writerow([Epoch1_ValidTrackNum, Epoch2_ValidTrackNum, Common_ValidTrackNum])
    fwriter.writerow([Epoch1_ERValidTrackNum, Epoch2_ERValidTrackNum, Common_ERValidTrackNum])
    fwriter.writerow([Epoch1_ValidTrackPortion, Epoch2_ValidTrackPortion, Common_ValidTrackPortion])
    fwriter.writerow([Epoch1_ERValidTrackPortion, Epoch2_ERValidTrackPortion, Common_ERValidTrackPortion])
    File.close()

    print('[Script][TimeCost]    :', time.perf_counter() - starttime0)
