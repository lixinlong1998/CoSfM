from config import *
import numpy as np
from src_CFTM import Func_CommonTiePoints
from src_Metashape import FuncMs_ErrorReduction

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


def loadFolder_CTPs(CTP_path):
    # read file and convert data format
    for root, dirs, files in os.walk(CTP_path):
        Matches = {}
        TracksReport = {}
        FilterTracksReport = {}
        CommonTracksReport = {}
        for fileName in files:
            filetype = fileName.split('_')[0]
            grid_id = eval(fileName[0:-4].split('_')[1])
            if not fileName[-3:] == 'npy':
                continue
            # open npy which not empty
            try:
                data_list = np.load(root + '/' + fileName, allow_pickle=True).tolist()
            except EOFError:
                continue
            # print(fileName)
            if filetype == 'Matches':
                Matches[grid_id] = data_list
            elif filetype == 'TracksReport':
                TracksReport[grid_id] = data_list
            elif filetype == 'FilterTracksReport':
                FilterTracksReport[grid_id] = data_list
            elif filetype == 'CommonTracksReport':
                CommonTracksReport[grid_id] = data_list
    return Matches, TracksReport, FilterTracksReport, CommonTracksReport


def readInnerMatches(MatchesReport):
    '''
    output:
        MatchesReport_e1 = {pair_id:[Matches, inlierMatches]}
            ,where Matches = [[keypoint_id, keypoint_id],...]
            ,where inlierMatches = [[keypoint_id, keypoint_id],...]
    '''

    def checkPairId(pair_id):
        if pair_id[0] > pair_id[1]:
            NewPair_id = (pair_id[1], pair_id[0])
        else:
            NewPair_id = pair_id
        return NewPair_id

    MatchesReport_e1 = {}
    MatchesReport_e2 = {}
    for grid_id, MatchesReport_cube in MatchesReport.items():
        MatchesReport_cube_e1, MatchesReport_cube_e2 = MatchesReport_cube[0], MatchesReport_cube[1]
        for MatchesReport_pair in MatchesReport_cube_e1:
            pair_id = checkPairId(MatchesReport_pair[0])  # 检查pair_id，小的数字在前,
            Matches = [[row[1], row[0]] for row in MatchesReport_pair[1]]
            InlierMatches = [[row[1], row[0]] for row in MatchesReport_pair[2]]
            # 按照pair_id 存储matches
            if pair_id not in MatchesReport_e1:
                MatchesReport_e1[pair_id] = [Matches, InlierMatches]  # Matches, inlierMatches
            else:
                MatchesReport_e1[pair_id][0] += Matches
                MatchesReport_e1[pair_id][1] += InlierMatches
        for MatchesReport_pair in MatchesReport_cube_e2:
            pair_id = checkPairId(MatchesReport_pair[0])  # 检查pair_id，小的数字在前
            Matches = [list(row) for row in MatchesReport_pair[1]]
            InlierMatches = [list(row) for row in MatchesReport_pair[2]]
            # 按照pair_id 存储matches
            if pair_id not in MatchesReport_e2:
                MatchesReport_e2[pair_id] = [Matches, InlierMatches]  # Matches, inlierMatches
            else:
                MatchesReport_e2[pair_id][0] += Matches
                MatchesReport_e2[pair_id][1] += InlierMatches

    return MatchesReport_e1, MatchesReport_e2


def reshapeMatchesStructure(MatchesReport_e1, MatchesReport_e2):
    '''
    output:
        MatchesSta = [InnerMatches1,InnerMatches2]
            ,where InnerMatches1 = [pair_MatchesNumber, ...]
        InlierMatchesSta = [InnerMatches1,InnerMatches2]
            ,where InnerMatches1 = [pair_MatchesNumber, ...]
    '''
    MatchesSta = [[], []]
    InlierMatchesSta = [[], []]
    for pair_id, MatchesReport_pair in MatchesReport_e1.items():
        Matches = MatchesReport_pair[0]
        inlierMatches = MatchesReport_pair[1]
        if len(Matches) > 0:
            MatchesSta[0].append(len(Matches))
        if len(inlierMatches) > 0:
            InlierMatchesSta[0].append(len(inlierMatches))
    for pair_id, MatchesReport_pair in MatchesReport_e2.items():
        Matches = MatchesReport_pair[0]
        inlierMatches = MatchesReport_pair[1]
        if len(Matches) > 0:
            MatchesSta[1].append(len(Matches))
        if len(inlierMatches) > 0:
            InlierMatchesSta[1].append(len(inlierMatches))
    return MatchesSta, InlierMatchesSta


def readCrossMatches(CommonTracksReport):
    '''
    output:
        CrossMatches = {pair_id:Matches}
            ,where Matches = [Match,...]
            ,where Match = [keypoint_id, keypoint_id,MatchQuality]
    '''
    CrossMatchesReport_dict = {}
    # 遍历所有的网格
    for grid_id, CommonTrackMatches_cube in CommonTracksReport.items():
        # 每个网格中有多个common track
        for CrossTrackMatches in CommonTrackMatches_cube:
            # 每个common track中对应了多个cross match
            for CrossMatch in CrossTrackMatches:
                feature1, feature2, MatchQuality = CrossMatch[0], CrossMatch[1], CrossMatch[2]
                if feature1[0] < feature2[0]:  # 小的数字在前
                    pair_id = (feature1[0], feature2[0])
                else:
                    pair_id = (feature2[0], feature1[0])
                Match = [feature1[1], feature2[1], MatchQuality]
                if pair_id not in CrossMatchesReport_dict:
                    CrossMatchesReport_dict[pair_id] = [Match]
                else:
                    CrossMatchesReport_dict[pair_id].append(Match)
    return CrossMatchesReport_dict


def readCommonTracksCrossMatch(CommonTracksReport):
    '''
    output:
        CTracksCrMatchReport = [CommonTrack_CrossMatches,...]
            ,where CommonTrack_CrossMatches = len(CrossTrackMatches)
    '''
    CTracksCrMatchReport = []
    for grid_id, CommonTrackReport in CommonTracksReport.items():
        # 每个网格中有多个common track
        for CrossTrackMatches in CommonTrackReport:
            CTracksCrMatchReport.append(len(CrossTrackMatches))
    return CTracksCrMatchReport


def readInnerTracks(TracksReport):
    '''
    output:
        InnerTracksReport_e1 = [InnerTrackReport,...]
            ,where InnerTrackReport = [viewsNum,MatchNum]
    '''
    InnerTracksReport_e1 = []
    InnerTracksReport_e2 = []
    for grid_id, TracksReport_cube in TracksReport.items():
        TracksReport_cube_e1, TracksReport_cube_e2 = TracksReport_cube[0], TracksReport_cube[1]
        for TracksReport in TracksReport_cube_e1:
            viewsNum = len(TracksReport[0])
            matchNum = len(TracksReport[1])
            InnerTracksReport_e1.append([viewsNum, matchNum])
        for TracksReport in TracksReport_cube_e2:
            viewsNum = len(TracksReport[0])
            matchNum = len(TracksReport[1])
            InnerTracksReport_e2.append([viewsNum, matchNum])
    return InnerTracksReport_e1, InnerTracksReport_e2


def analyseTracks(TracksEpoch, points, point_ids):
    TracksReport = []
    for track_id, epochNum in enumerate(TracksEpoch):
        # TrackType
        signal = Func_CommonTiePoints.analyseSignal_Epoch(epochNum)
        if sum(signal) > 1:
            TrackType = -1
        else:
            for epoch_id, sig in enumerate(signal):
                if sig > 0:
                    TrackType = epoch_id  # 0,1,2,3...
        # point_id
        point_id = point_ids[track_id]
        # point_valid
        if point_id == -1:
            point_valid = 0
        else:
            if points[point_id].valid:
                point_valid = 1
            else:
                point_valid = 0
        # track_id, point_id, TrackType, point_valid, point_ERvalid(set all points valid before ER as default)
        TracksReport.append([track_id, point_id, TrackType, point_valid, point_valid])
    return TracksReport


def analyseRemovedTracks(TracksReport, points, removedPoints):
    # assign 0 to the removed points after ER
    for point_id in removedPoints:
        track_id = points[point_id].track_id
        TracksReport[track_id][4] = 0
    # split TrackReport by epoch category
    Epoch1_TrackReport = []
    Epoch2_TrackReport = []
    Common_TrackReport = []
    TrackTypes = sorted(set([TrackReport[2] for TrackReport in TracksReport]))  # from -1 to 0,1,2,3,...
    for track_id, TrackReport in enumerate(TracksReport):
        if TrackReport[2] == TrackTypes[0]:
            Common_TrackReport.append(TrackReport)
        if TrackReport[2] == TrackTypes[1]:
            Epoch1_TrackReport.append(TrackReport)
        if TrackReport[2] == TrackTypes[2]:
            Epoch2_TrackReport.append(TrackReport)
    return Epoch1_TrackReport, Epoch2_TrackReport, Common_TrackReport


def getTracksReport(chunk, TracksEpoch, points, point_ids, CriteriasThreshold=[0.5, 20, 8]):
    '''
    input:
        CriteriasThreshold = [0.3, 15, 5]  # RE,RU,PA
    output:
         Epoch1_TrackReport = [TrackReport,...]
         Epoch2_TrackReport = [TrackReport,...]
         Common_TrackReport = [TrackReport,...]
            ,where TrackReport = [track_id, point_id, TrackType, point_valid, point_ERvalid]
    '''
    # analysis track
    TracksReport = analyseTracks(TracksEpoch, points, point_ids)
    # Error Reduction but don't update, the valid attribute of removed points  will be assigned a False
    chunk, removedPoints = FuncMs_ErrorReduction.ErrorReduction(chunk, CriteriasThreshold, update=False)
    # analysis removed Points
    return analyseRemovedTracks(TracksReport, points, removedPoints)


def getMatchAnalyse(chunkMatches, chunkDescriptors, camera_ids, CamerasEpoch, epochs):
    '''
    input:
    output:
        MatchesSta = [[n,n,...,n],[n,n,...,n],[n,n,...,n]]
            ,with three list contains the pairs in epoch1,epoch2 and cross-epoch respectively
            ,where n is matches number of pair
        MatchesDis = [[pair,...],[pair,...],[pair,...]]
            ,where pair =  [distance,distance,...,distance]
    '''
    MatchesSta = [[], [], []]  # epoch1,epoch2,cross-epoch
    MatchesDis = [[], [], []]  # epoch1,epoch2,cross-epoch
    # 根据epochs的数量创建空列表
    EpochMatchesSta = [[] for i in range(len(epochs))]
    EpochMatchesDis = [[] for i in range(len(epochs))]
    for pair_id, matches in chunkMatches.items():
        identifier1, identifier2 = str(pair_id[0]), str(pair_id[1])
        # 查看图像属于哪个时相
        image1_epoch = CamerasEpoch[camera_ids[identifier1]]
        image2_epoch = CamerasEpoch[camera_ids[identifier2]]
        # 计算每个match的distance
        Descriptors1 = chunkDescriptors[int(identifier1)][matches[:, 0]]
        Descriptors2 = chunkDescriptors[int(identifier2)][matches[:, 1]]
        distances = np.linalg.norm(Descriptors1 - Descriptors2, axis=1)
        # 检查该像对是否跨时相
        if image1_epoch != image2_epoch:
            MatchesSta[2].append(len(matches))
            MatchesDis[2].append(list(distances))
        else:
            EpochMatchesSta[image1_epoch].append(len(matches))
            EpochMatchesDis[image1_epoch].append(list(distances))
    # 将innerEpochMatches中不为空的时相匹配数列表依次传递给MatchesSta的前两个空列表
    EpochMatches_notnone = [MatchesNumList for MatchesNumList in EpochMatchesSta if MatchesNumList]
    MatchesSta[0] = EpochMatches_notnone[0]
    MatchesSta[1] = EpochMatches_notnone[1]
    # 将innerEpochMatches中不为空的时相匹配数列表依次传递给MatchesDis的前两个空列表
    EpochMatches_notnone = [MatchesNumList for MatchesNumList in EpochMatchesDis if MatchesNumList]
    MatchesDis[0] = EpochMatches_notnone[0]
    MatchesDis[1] = EpochMatches_notnone[1]
    return MatchesSta, MatchesDis


def alignAllMatch(chunkMatches, chunkInlierMatches, camera_ids, CamerasEpoch, epochs):
    """
    output:
        MatchesReport_Epoch1 = {Pair:[MatchesNum, inlierMatchesNum],...}
        MatchesReport_Epoch2 = {Pair:[MatchesNum, inlierMatchesNum],...}
        MatchesReport_Common = {Pair:[MatchesNum, inlierMatchesNum],...}
    """
    MatchesReport_Common = {}
    MatchesReport_inner = [{} for i in range(len(epochs))]
    for pair_id, matches in chunkMatches.items():
        identifier1, identifier2 = str(pair_id[0]), str(pair_id[1])
        # 查看图像属于哪个时相
        image1_epoch = CamerasEpoch[camera_ids[identifier1]]
        image2_epoch = CamerasEpoch[camera_ids[identifier2]]
        # 检查该像对是否跨时相
        if image1_epoch != image2_epoch:
            MatchesReport_Common[pair_id] = [len(matches), 0]
        else:
            MatchesReport_inner[image1_epoch][pair_id] = [len(matches), 0]
    # 将innerEpochMatches中不为空的时相匹配数列表依次传递给MatchesSta的前两个空列表
    MatchesReport_inner = [MatchesReport for MatchesReport in MatchesReport_inner if len(MatchesReport) != 0]
    MatchesReport_Epoch1 = MatchesReport_inner[0]
    MatchesReport_Epoch2 = MatchesReport_inner[1]
    # 把inlier的信息也写进来
    for pair_id, InlierMatches in chunkInlierMatches.items():
        if pair_id in MatchesReport_Common:
            MatchesReport_Common[pair_id][1] = len(InlierMatches)
        elif pair_id in MatchesReport_Epoch1:
            MatchesReport_Epoch1[pair_id][1] = len(InlierMatches)
        elif pair_id in MatchesReport_Epoch2:
            MatchesReport_Epoch2[pair_id][1] = len(InlierMatches)
    return MatchesReport_Epoch1, MatchesReport_Epoch2, MatchesReport_Common


def getPairInlierMatchRatio(MatchesReport):
    """
    input:
        MatchesReport = {Pair:[MatchesNum, inlierMatchesNum],...}
    output:
        MatchesReport = {Pair:[MatchesNum, inlierMatchesNum],...}
    """
    PIMR = []
    for pair_id, AllMatchesNum in MatchesReport.items():
        # 跳过无内点匹配的像对
        if AllMatchesNum[1] == 0:
            continue
        PairInlierMatchRatio = AllMatchesNum[1] / AllMatchesNum[0]
        PIMR.append(PairInlierMatchRatio)
    return np.mean(PIMR)
