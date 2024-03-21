import networkx as nx
import numpy as np

import random
import csv
from scipy.special import comb
import src_CFTM.Func_Triangulate as Tri
import src_CFTM.Func_Camera as Cam
import src_CFTM.Func_Match as Matc


#################################################   generateTracks  ####################################################

def generateTracks(matches_cube):
    '''
    input:
        matches_cube = [matches_cube_pair_id,matches_cube_pair_id,...,matches_cube_pair_id]
            ,where matches_cube_pair_id = [(camera_id1, camera_id2),match,match,...,match]
            ,where match = [[camera_id, projection_info, keypoints_id],[camera_id, projection_info, keypoints_id], matchQuality]
            ,where projection_info = [u,v]
            ,where matchQuality = [distance, point3D]
    output:
        Tracks_Generated = [Track,Track,...,Track]
            ,where Track = set{feature,feature,...,feature}
            ,where feature = tuple(camera_id, keypoints_id)
        Tracks_Matches = [Track_Matches,Track_Matches,...,Track_Matches]
            ,where Track_Matches = [Match,Match,...,Match]
            ,where Match = (feature1, feature2, 1)
            ,where feature = tuple(camera_id, keypoints_id)
    '''
    nodes = getNodes(matches_cube)
    edges = getEdges(matches_cube)
    Tracks_Generated, Tracks_Matches = ConnectedComponent(nodes, edges)
    return Tracks_Generated, Tracks_Matches


def getNodes(matches_cube):
    '''
    input:
        matches_cube = [matches_cube_pair_id,matches_cube_pair_id,...,matches_cube_pair_id]
            ,where matches_cube_pair_id = [(camera_id1, camera_id2),match,match,...,match]
            ,where match = [[camera_id, projection_info, keypoints_id],[camera_id, projection_info, keypoints_id], matchQuality]
            ,where projection_info = [u,v]
            ,where matchQuality = [distance, point3D]
    output:
        nodes = [(camera_id, keypoints_id),...,(camera_id, keypoints_id)]
    note: ConnectedComponent didn't accept the projection_info in the nodes, so we don't use it here.
    '''
    # The matched keypoints are defined as nodes
    nodes = []
    for matches_cube_pair_id in matches_cube:
        camera_id1 = matches_cube_pair_id[0][0]
        camera_id2 = matches_cube_pair_id[0][1]
        for match in matches_cube_pair_id[1:]:
            nodes.append((camera_id1, match[0][2]))
            nodes.append((camera_id2, match[1][2]))
    return list(set(nodes))


def getEdges(matches_cube):
    '''
    input:
        matches_cube = [matches_cube_pair_id,matches_cube_pair_id,...,matches_cube_pair_id]
            ,where matches_cube_pair_id = [(camera_id1, camera_id2),match,match,...,match]
            ,where match = [[camera_id, projection_info, keypoints_id],[camera_id, projection_info, keypoints_id], matchQuality]
            ,where projection_info = [u,v]
            ,where matchQuality = [distance, point3D]
    output:
        edges = ((camera_id1, keypoints_id1), (camera_id2, keypoints_id2), 1),...)
    '''
    edges = []
    for matches_cube_pair_id in matches_cube:
        camera_id1 = matches_cube_pair_id[0][0]
        camera_id2 = matches_cube_pair_id[0][1]
        for match in matches_cube_pair_id[1:]:
            # edge is a tuple with element1,2 and weight (here is 1)
            edges.append(((camera_id1, match[0][2]), (camera_id2, match[1][2]), 1))
    return edges


def ConnectedComponent(nodes, edges):
    '''
    input:
        nodes = [feature,feature,...,feature]
        edges = [(feature1, feature2, 1),...]
            ,where feature = tuple(camera_id, keypoints_id)
    output:
        ConnectedComponent = [Track,Track,...,Track]
            ,where Track = set{feature,feature,...,feature}
            ,where feature = tuple(camera_id, keypoints_id)
        SubGraphsEdge = [subGraph_edges,subGraph_edges,...,subGraph_edges]
            ,where subGraph_edges = [edge,edge,...,edge]
            ,where edge = (feature1, feature2, 1)
            ,where feature = tuple(camera_id, keypoints_id)
    '''
    # generate connected component
    Graph = nx.Graph()
    Graph.add_nodes_from(nodes)
    Graph.add_weighted_edges_from(edges)
    ConnectedComponent = sorted(nx.connected_components(Graph), key=len, reverse=True)

    # get matches of track
    SubGraphsEdge = []
    for subGraph_nodes in ConnectedComponent:
        subGraph = Graph.subgraph(subGraph_nodes)
        SubGraphsEdge.append(subGraph.edges())
    return ConnectedComponent, SubGraphsEdge


###################################################  filterTracks  #####################################################

def filterTracks(grid_id, Tracks_Generated, Cameras, Sensors, chunkKeypoints, ratio_inlier_init, confidence,
                 criterias, threshold_RepError2, show=True):
    '''RANSAC tracks filter
    input:
        Tracks_Generated = [TrackG,TrackG,...,TrackG]
            ,where TrackG = set{feature,feature,...,feature}
            ,where feature = tuple(camera_id, keypoints_id)
        Cameras =
        Sensors =
        ratio_inlier_init = An initial guess of the average ratio of inliers in track, here set as 0.03
        confidence = A probability of being able to get the consensus set by randomly sample certain times
        criterias = [threshold_RepError, threshold_Angle, *threshold_ProAcc], here set as [8, 2]
        threshold_RepError2 = A threshold used for the final assessment of point quality of track
    output:
        Tracks_Filtered = [Track_Filtered,Track_Filtered,...,Track_Filtered]
            ,where Track_Filtered = [view,view,...,view]
            ,where view = [camera_id, projection_info, index_info]
            ,where projection_info = [u, v]
            ,where index_info = [identifier, keypoints_id, originTrack_id]
        Tracks_Quality = [Track_Quality,Track_Quality,...,Track_Quality]
            ,where Track_Quality = [Point3D, PointRepError]
            ,where point3D = ndarray([X,Y,Z,1])
            ,where PointRepError = [ReprojectError, ReprojectError,..., ReprojectError]

    algorithm:
        这个函数可以使用ransac、三角交会、匹配数来过滤tracks，目的是尽可能保留具有良好交会条件和时相匹配条件的tracks，并记录它们的时相匹配条件，以用作后续的排序和优选。
        step1  通过读取相机的epoch，将track的view_ids分为两个子集：Track_ViewId_epoch1 &Track_ViewId_epoch2
        step2  基于子集中view_id的数量，过滤掉那些view_id少于2的track
        step3
    '''
    # [test part]
    testFile = open('G:/AResearchG/20221223_CoSfM/Code/test/Test_RANSACtrackfiltter/test_filterTracks.txt', "w")
    fwriter = csv.writer(testFile, delimiter='\t', lineterminator='\n')
    Tracks_Generated_Length_Statistic = []

    # parameters
    threshold_TrackLength = 2
    threshold_Eligible = 3
    ratio_inlier_AVG = ratio_inlier_init
    ratio_inlier_List = []

    # output list
    removeTrackIdList = []
    Tracks_Filtered = []
    Tracks_Quality = []

    for Track_id, TrackG in enumerate(Tracks_Generated):
        Tracks_Generated_Length_Statistic.append(len(TrackG))
        # 1  filter: Only keep the Track length bigger than threshold
        #              this condition guarantee that there is at least one match across epochs
        if len(TrackG) < threshold_TrackLength:
            removeTrackIdList.append(Track_id)
            if show:
                print('[Script]    [{}]    remove track:{}, too less views in track'.format(grid_id, Track_id))
            continue
        else:
            Track = convertTrackData(TrackG, chunkKeypoints, Cameras, Track_id)

        # 2  filter: Using RANSAC to filter track
        if len(Track) >= 3:
            Track_filtered_List, ratio_inliers = filterTrack_RANSAC(Track, Cameras, Sensors, ratio_inlier_AVG,
                                                                    threshold_Eligible, confidence, criterias)
            # Check the results of RANSAC
            if len(Track_filtered_List) == 0:
                # It means that there is no consensus set in this track
                if show:
                    print('[Script]    [{}]    remove track:{}, no consensus could be found'.format(grid_id, Track_id))
                continue
            else:
                # The inlier ratio corresponding to the consensus sets is added to the external statistical list.
                ratio_inlier_List += ratio_inliers
                ratio_inlier_AVG = np.mean(ratio_inlier_List)
                # When the track is very short, the consensus is usually the track itself, resulting in an inlier rate of 1
                # We try to avoid this case that may leads to very few iterations
                if ratio_inlier_AVG > 0.5:
                    ratio_inlier_AVG = 0.5
            for Track in Track_filtered_List:
                # 3  filter: Using points quality to filter track
                Point3D, PointRepError = filterTrack_Point(Track, Cameras, Sensors)
                if Point3D[3] == 0:
                    # print('[Script]    [{}]    remove track:{}, Triangulate failed'.format(grid_id, Track_id))
                    continue
                elif max(PointRepError) > threshold_RepError2:
                    # print('[Script]    [{}]    remove track:{}, Reprojection error too large'.format(grid_id, Track_id))
                    continue
                # 4  Output filtered track
                Tracks_Filtered.append(Track)
                Tracks_Quality.append([Point3D, PointRepError])
        else:
            # 3  filter: Using points quality to filter track
            Point3D, PointRepError = filterTrack_Point(Track, Cameras, Sensors)
            if Point3D[3] == 0:
                removeTrackIdList.append(Track_id)
                if show:
                    print('[Script]    [{}]    remove track:{}, Triangulate failed'.format(grid_id, Track_id))
                continue
            elif max(PointRepError) > threshold_RepError2:
                removeTrackIdList.append(Track_id)
                if show:
                    print('[Script]    [{}]    remove track:{}, Reprojection error too large'.format(grid_id, Track_id))
                continue
            # 4  Output filtered track
            Tracks_Filtered.append(Track)
            Tracks_Quality.append([Point3D, PointRepError])

    print('max track length:', max(Tracks_Generated_Length_Statistic))
    print('average track length:', np.mean(Tracks_Generated_Length_Statistic))
    print('min track length:', min(Tracks_Generated_Length_Statistic))
    print('remove tracks:', len(removeTrackIdList))
    print('remain tracks:', len(Tracks_Filtered))
    return Tracks_Filtered, Tracks_Quality


def convertTrackData(TrackG, chunkKeypoints, Cameras, originTrack_id):
    '''
    input:
        TrackG = set{feature,feature,...,feature}
            ,where feature = tuple(camera_id, keypoints_id)
    output:
        Track = [view,view,...view]
            ,where view = [camera_id, projection_info, index_info]
            ,where projection_info = [u, v]
            ,where index_info = [identifier, keypoints_id, originTrack_id]
    '''
    Track = []
    for camera_id, keypoints_id in list(TrackG):
        identifier = Cameras[camera_id][-1]
        projection_info = list(chunkKeypoints[identifier][keypoints_id, 0:2])
        view = [camera_id, projection_info, [identifier, keypoints_id, originTrack_id]]
        Track.append(view)
    return Track


def adaptiveIteration(ratio, confidence):
    '''
    input:
        length = track_Length_AVG
        ratio = ratio_inlier
        confidence = probability in (0,1), usually set as 0.95 or 0.99
    output:
        number of iterations
    algorithm:
        Iterations = ceil[ln(1 - confidence) / ln(1 - ratio^2)]
    '''
    iterations = np.log(1 - confidence) / np.log(1 - ratio * ratio)
    return int(iterations) + 1


def RANSAC(track, num_iterations, Cameras, Sensors, threshold_Eligible, threshold_Angle, threshold_RepError1):
    def checkFromSameCamera(ViewId_sample):
        camera_id1 = track[ViewId_sample[0]][0]
        camera_id2 = track[ViewId_sample[1]][0]
        if camera_id1 == camera_id2:
            return True
        else:
            return False

    assert len(track) >= 3
    best_num_inliers = 2
    max_sample_combinations = comb(len(track), 2)
    ViewId_sample_list = []
    best_ViewId_inlier = []
    best_ViewId_outlier = []
    for i in range(num_iterations):
        # Randomly select a subset of correspondences
        while True:
            ViewId_sample = random.sample(range(len(track)), 2)
            if set(ViewId_sample) in ViewId_sample_list:
                if len(ViewId_sample_list) == max_sample_combinations:
                    return [track[s] for s in best_ViewId_inlier], [track[t] for t in best_ViewId_outlier]
                else:
                    continue
            else:
                if checkFromSameCamera(ViewId_sample):
                    ViewId_sample_list.append(set(ViewId_sample))
                    continue
                else:
                    ViewId_sample_list.append(set(ViewId_sample))
                    break

        # Estimate the 3D positions of the corresponding points using linear triangulation
        track_sample = [track[j] for j in ViewId_sample]
        Point3D_sample = Tri.Triangulate(track_sample, Cameras, Sensors)
        # Check1: Triangulate successfully?
        if Point3D_sample[3] == 0:
            # print('[Script]        Triangulate failed!')
            continue
        # Check2: Triangulation angle big enough?
        angle = Tri.TriangulationAngle(Cameras[track_sample[0][0]], Cameras[track_sample[1][0]], Point3D_sample,
                                       dtype='deg')
        if angle < threshold_Angle:
            # print('[Script]        Failed: the triangulation angle is too small!')
            continue

        # Compute the reprojection errors for all views in track
        Point3D_sample_RepError = np.zeros(len(track))
        for view_id, view in enumerate(track):
            Point3D_sample_RepError[view_id] = Cam.ReprojectError(
                Point3D_sample[:-1], Cameras[view[0]], Sensors, [view[1][0], view[1][1]], 'L1')
        # Count the number of inliers (correspondences with reprojection errors below the threshold)
        ViewId_inlier = np.where(Point3D_sample_RepError < threshold_RepError1)[0]  # array([])
        ViewId_outlier = np.where(Point3D_sample_RepError >= threshold_RepError1)[0]  # array([])

        if len(ViewId_inlier) > best_num_inliers:
            # Update the best subset of inliers found so far
            best_num_inliers = len(ViewId_inlier)
            best_ViewId_inlier = ViewId_inlier
            best_ViewId_outlier = ViewId_outlier
            # Check if this consensus is already satisfy the inlier ratio threshold
            if best_num_inliers >= threshold_Eligible or len(ViewId_outlier) == 0:
                return [track[s] for s in best_ViewId_inlier], [track[t] for t in best_ViewId_outlier]
            else:
                continue
        else:
            continue
    return [track[s] for s in best_ViewId_inlier], [track[t] for t in best_ViewId_outlier]


def filterTrack_RANSAC(Track, Cameras, Sensors, ratio_inlier_AVG, threshold_Eligible, confidence, criterias):
    '''
    input:
        Track = [view,view,...view]
            ,where view = [camera_id, projection_info, index_info]
            ,where projection_info = [u, v]
            ,where index_info = [identifier, keypoints_id, originTrack_id]
        Cameras
        Sensors
        ratio_inlier_AVG
        threshold_Eligible
        confidence = 0.9
        criterias = [threshold_RepError, threshold_Angle, *threshold_ProAcc]
    output:
        Track_filtered_List = [Track,Track,...,Track]
            ,where Track = [view,view,...view]
            ,where view = [camera_id, projection_info, index_info]
            ,where projection_info = [u, v]
            ,where index_info = [identifier, keypoints_id, originTrack_id]
        ratio_inlier_List = [ratio,...]
    algorithm:
    '''
    # get value
    threshold_RepError1 = criterias[0]
    threshold_Angle = criterias[1]
    Track_Length_initial = len(Track)

    # Recursively RANSAC
    ratio_inlier_List = []
    Track_filtered_List = []
    while True:
        # Use adaptive iteration number
        num_iterations = adaptiveIteration(ratio_inlier_AVG, confidence)
        # Implement RANSAC
        Track_Consensus, Track = RANSAC(Track, num_iterations, Cameras, Sensors,
                                        threshold_Eligible, threshold_Angle, threshold_RepError1)
        if len(Track_Consensus) > 0:
            # When a consensus set exists, the inlier ratio corresponding to this consensus set is recorded
            # But when the track is small, consensus is likely to be the original track, where ratio_inlier is 1
            # We try to avoid this case that may leads to very few iterations
            ratio_inlier_List.append(len(Track_Consensus) / Track_Length_initial)
            Track_filtered_List.append(Track_Consensus)
            if len(Track_Consensus) > 2 and len(Track) > 2:
                # try to find the consensus set from the remaining measurements
                continue
            else:
                return Track_filtered_List, ratio_inlier_List
        else:
            # no consensus
            return [], []
    return Track_filtered_List, ratio_inlier_List


def filterTrack_Point(Track, Cameras, Sensors):
    '''
    input:
        Track = [view,view,...view]
            ,where view = [camera_id, projection_info, index_info]
            ,where projection_info = [u, v]
            ,where index_info = [identifier, keypoints_id, originTrack_id]
        Cameras
        Sensors
    output:
        Point_epoch = ndarray([X,Y,Z,1])
        PointRepError_epoch = [ReprojectError, ReprojectError,..., ReprojectError]
    algorithm:
        三角交会后得到三维点，并重投影回可视像片，计算并统计重投影误差
    '''
    Point3D = Tri.Triangulate(Track, Cameras, Sensors)
    # Compute the reprojection errors for all correspondences
    PointRepError_epoch = []
    for view_id, view in enumerate(Track):
        PointRepError_epoch.append(
            Cam.ReprojectError(Point3D[:-1], Cameras[view[0]], Sensors, [view[1][0], view[1][1]], 'L1'))
    return Point3D, PointRepError_epoch


############################################  Match common tracks  #####################################################
def matchTracks(Tracks_epoch1, Tracks_epoch2, Tracks_Quality_epoch1, Tracks_Quality_epoch2, Cameras, Sensors,
                chunkDescriptors, threshold_Distance, crossCheck=True):
    '''match tracks across epochs
    input:
        Tracks_epoch = [Track,Track,...,Track]
            ,where Track = [view,view,...,view]
            ,where view = [camera_id, projection_info, index_info]
            ,where projection_info = [u, v]
            ,where index_info = [identifier, keypoints_id, originTrack_id]
        Tracks_Quality_epoch = [Track_Quality,Track_Quality,...,Track_Quality]
            ,where Track_Quality = [Point3D, PointRepError]
            ,where point3D = ndarray([X,Y,Z,1])
            ,where PointRepError = [ReprojectError, ReprojectError,..., ReprojectError]
        Cameras =
        chunkDescriptors =
        threshold_Distance = The matching range is limited by the coordinates of the 3D points
    output:
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
    algorithm:

    '''

    def findNearestNeighbours(similarity_table_track1id, Tracks_epoch1):
        '''
        output:
            NearestNeighbours_track1id=[track2id,-1,track2id,...] indexed by track1id,where -1 means no track2 near to it
        '''
        NearestNeighbours_track1id = [-1 for i in range(len(Tracks_epoch1))]
        for Track1_id, Track1_neighbours in enumerate(similarity_table_track1id):
            if not Track1_neighbours:
                continue
            elif len(Track1_neighbours) == 1:
                NearestNeighbours_track1id[Track1_id] = Track1_neighbours[0][0]
            else:
                TN = np.asarray(Track1_neighbours)
                TN_ranked = TN[np.argsort(TN[:, 1])]  # rank similarity from min to max
                NearestNeighbours_track1id[Track1_id] = int(TN_ranked[0, 0])  # Track2_id
        return NearestNeighbours_track1id

    CommonTracks = []
    CommonTracksMatches = []
    if crossCheck:
        CommonTracks_tentative = {}
        CommonTracksMatches_tentative = {}
        similarity_table_track1id = [[] for i in range(len(Tracks_epoch1))]
        similarity_table_track2id = [[] for i in range(len(Tracks_epoch2))]
        for Track1_id, Track1 in enumerate(Tracks_epoch1):
            Point3D1 = Tracks_Quality_epoch1[Track1_id][0]
            # PointRepError1 = Tracks_Quality_epoch1[Track1_id][-1]
            for Track2_id, Track2 in enumerate(Tracks_epoch2):
                Point3D2 = Tracks_Quality_epoch2[Track2_id][0]
                # PointRepError2 = Tracks_Quality_epoch2[Track2_id][-1]
                # The matching range is limited by the coordinates of the 3D points
                if np.linalg.norm(Point3D2[:-1] - Point3D1[:-1]) > threshold_Distance:
                    continue
                # Use brute force matcher to calculate the matches between track1 and track2.
                # this always matched successfully, due to the use of cross-check in outerMatch
                CommonTrackMatches = Matc.outerMatch(Track1, Track2, Cameras, Sensors, chunkDescriptors)
                # [Exception Handling]: In case there is no epoch match in Track pair.
                if len(CommonTrackMatches) == 0:
                    # skip this Track pair
                    continue
                # Use the minimum distance of cross-match to represent the similarity of two given tracks
                similarity = min([M[2][0] for M in CommonTrackMatches])
                similarity_table_track1id[Track1_id].append([Track2_id, similarity])
                similarity_table_track2id[Track2_id].append([Track1_id, similarity])
                # for matched Track pair, align them and append to tentative result
                CommonTracks_tentative[(Track1_id, Track2_id)] = Track1 + Track2
                CommonTracksMatches_tentative[(Track1_id, Track2_id)] = CommonTrackMatches
        # find the Nearest Neighbours for each track
        NearestNeighbours_track1id = findNearestNeighbours(similarity_table_track1id, Tracks_epoch1)
        NearestNeighbours_track2id = findNearestNeighbours(similarity_table_track2id, Tracks_epoch2)
        # cross check for getting goodmatches
        for Track1_id, Track2_id in enumerate(NearestNeighbours_track1id):
            if Track2_id == -1:
                continue
            elif NearestNeighbours_track2id[Track2_id] == Track1_id:
                # good matches
                CommonTracks.append(CommonTracks_tentative[(Track1_id, Track2_id)])
                CommonTracksMatches.append(CommonTracksMatches_tentative[(Track1_id, Track2_id)])
    else:
        for Track1_id, Track1 in enumerate(Tracks_epoch1):
            Point3D1 = Tracks_Quality_epoch1[Track1_id][0]
            test_matchingTrackNum = 0
            test_matchedTrackNum = 0
            for Track2_id, Track2 in enumerate(Tracks_epoch2):
                Point3D2 = Tracks_Quality_epoch2[Track2_id][0]
                # The matching range is limited by the coordinates of the 3D points
                if np.linalg.norm(Point3D2[:-1] - Point3D1[:-1]) > threshold_Distance:
                    continue
                test_matchingTrackNum += 1
                # Use brute force matcher to calculate the matches between track1 and track2
                # this always matched successfully, due to the use of cross-check in outerMatch
                CommonTrackMatches = Matc.outerMatch(Track1, Track2, Cameras, Sensors, chunkDescriptors)
                # [Exception Handling]: In case there is no epoch match in Track pair.
                if len(CommonTrackMatches) == 0:
                    # skip this Track pair
                    continue
                test_matchedTrackNum += 1
                # for matched Track pair, align them and append to result
                CommonTracks.append(Track1 + Track2)
                CommonTracksMatches.append(CommonTrackMatches)
    return CommonTracks, CommonTracksMatches


############################################ Ccommon tracks Filter #####################################################
def CommonTrackRatioFilter(CommonTracks, CommonTracksMatches, ratio=0.7):
    '''
    input/output:
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
    '''


################################################  Select tracks  #######################################################
# def SelectBestTracks(CommonTracks, CommonTracksReport, Numbers):
#     '''match tracks across epochs
#     input:
#         CommonTracks = [Track_Matched,Track_Matched,...,Track_Matched]
#             ,where Track_Matched = [view,view,...,view]
#             ,where view = [camera_id, projection_info, index_info]
#             ,where projection_info = [u, v]
#             ,where index_info = [identifier, keypoints_id, originTrack_id]
#         CommonTracksReport = [Track_Report,Track_Report,...,Track_Report]
#             ,where Track_Report = [similarity, angles, RepErrors]
#             ,where similarity = [distance,distance,...,distance],distance is the L2 normal distance between descriptors of match.
#             ,where angles = [angle,angle,...,angle], angle is the triangulation angle between epoch match.
#             ,where RepErrors = [RepError, RepError,..., RepError], RepError is calculated separately from epochs
#         Numbers = [N1, N2, N3, N4]
#             ,where N1 = 匹配数是跨时相特征点控制配准效果的基本参数
#             ,where N2 = 用最大距离来尽可能过滤掉包含误匹配的track
#             ,where N3 = 有较好匹配的点很有可能在空间中很接近，导致过小的交会角，这不利于像片位姿的对齐
#             ,where N4 =
#     output:
#         TracksMatched_Best = [Track,Track,...,Track]
#             ,where Track = [view,view,...,view]
#             ,where view = [camera_id, projection_info, keypoints_id, distance]
#             ,where projection_info = [u, v]
#         TracksReport_Best = [Track_Report,Track_Report,...,Track_Report]
#             ,where Track_Report = [similarity, angles, RepErrors]
#             ,where similarity = [distance,distance,...,distance],distance is the L2 normal distance between descriptors of match.
#             ,where angles = [angle,angle,...,angle], angle is the triangulation angle between epoch match.
#             ,where RepErrors = [RepError, RepError,..., RepError], RepError is calculated separately from epochs
#     algorithm:
#     '''
#     # construct Tracks_Quality rank list
#     Tracks_Quality_index = []
#     for Track_id, Track_Report in enumerate(CommonTracksReport):
#         # read report
#         similarity = Track_Report[0]
#         angles = Track_Report[1]
#         RepErrors = Track_Report[2]
#         # indicators
#         Num_Matches = len(similarity)
#         Distance_Max = max(similarity)
#         Angles_Min = min(angles)
#         RepError_Max = max(RepErrors)
#         Tracks_Quality_index.append([Track_id, Num_Matches, Distance_Max, Angles_Min, RepError_Max])
#     Indexes = np.asarray(Tracks_Quality_index)
#
#     def getFirstNRows(Indexes_ranked, N):
#         if len(Indexes_ranked) > N:
#             Indexes_sorted = Indexes_ranked[:N, :]
#         else:
#             Indexes_sorted = Indexes_ranked
#         return Indexes_sorted
#
#     ''' select top N1 tracks by epoch matches number
#     匹配数是跨时相特征点控制配准效果的基本参数
#     '''
#     N1 = Numbers[0]
#     Indexes_ranked = Indexes[np.argsort(-Indexes[:, 1])]  # from max to min
#     Indexes = getFirstNRows(Indexes_ranked, N1)
#     ''' select top N2 tracks by Distance_Max
#     光看匹配数不行，因为如果匹配数过多，很可能意味着存在误匹配的概率大
#     用最大距离来尽可能过滤掉包含误匹配的track
#     '''
#     N2 = Numbers[1]
#     Indexes_ranked = Indexes[np.argsort(Indexes[:, 2])]  # from min to max
#     Indexes = getFirstNRows(Indexes_ranked, N2)
#     ''' select top N3 tracks by across angle of view rays
#     有较好匹配的点很有可能在空间中很接近，导致过小的交会角，这不利于像片位姿的对齐
#     '''
#     N3 = Numbers[2]
#     Indexes_ranked = Indexes[np.argsort(-Indexes[:, 3])]  # from max to min
#     Indexes = getFirstNRows(Indexes_ranked, N3)
#     ''' select top N4 tracks by Distance_Max
#     除了考虑匹配效果，还需要通过重投影误差来判断特征点是否位于稳定区域，因为在非稳定区域（例如植被、河流）的点通常具有较大的重投影误差
#     '''
#     N4 = Numbers[3]
#     Indexes_ranked = Indexes[np.argsort(Indexes[:, 4])]  # from min to max
#     Indexes = getFirstNRows(Indexes_ranked, N4)
#     # output Tracks_Best
#     TracksMatched_Best = [CommonTracks[int(Index[0])] for Index in Indexes]
#     TracksReport_Best = [CommonTracksReport[int(Index[0])] for Index in Indexes]
#     return TracksMatched_Best, TracksReport_Best


################################################  Report tracks  #######################################################
def getOriginTracksReport(Tracks_Generated, Tracks_Matches):
    '''
    input:
        Tracks_Generated = [Track,Track,...,Track]
            ,where Track = set{feature,feature,...,feature}
            ,where feature = tuple(camera_id, keypoints_id)
        Tracks_Matches = [Track_Matches,Track_Matches,...,Track_Matches]
            ,where Track_Matches = [Match,Match,...,Match]
            ,where Match = (feature1, feature2, 1)
            ,where feature = tuple(camera_id, keypoints_id)
    output:
        TracksReport_cube_e1=[TrackReport, TrackReport,…, TrackReport]
            ,where TrackReport = [views,edges]
            ,where views = np.array([[camera_id, keypoint_id],…])
            ,where edges = np.array([[view_id, view_id],…])
    '''
    TracksReport_cube = []
    for Track_id, Track in enumerate(Tracks_Generated):
        Track = list(Track)
        edges = []
        for Match in Tracks_Matches[Track_id]:
            view1_id = Track.index(Match[0])
            view2_id = Track.index(Match[1])
            edges.append([view1_id, view2_id])
        views = [[feature[0], feature[1]] for feature in Track]
        TracksReport_cube.append([np.asarray(views), np.asarray(edges)])
    return TracksReport_cube


def getFilterTracksReport(Tracks_Filtered, Tracks_Matches):
    '''
    input:
        Tracks_Filtered = [Track_Filtered,Track_Filtered,...,Track_Filtered]
            ,where Track_Filtered = [view,view,...,view]
            ,where view = [camera_id, projection_info, index_info]
            ,where projection_info = [u, v]
            ,where index_info = [identifier, keypoints_id, originTrack_id]
        Tracks_Matches = [Track_Matches,Track_Matches,...,Track_Matches]
            ,where Track_Matches = [Match,Match,...,Match]
            ,where Match = (feature1, feature2, 1)
            ,where feature = tuple(camera_id, keypoints_id)
    output:
        TracksReport_cube_e1=[TrackReport, TrackReport,…, TrackReport]
            ,where TrackReport = [views,edges]
            ,where views = np.array([[camera_id, keypoint_id],…])
            ,where edges = np.array([[view_id, view_id],…])
    '''
    TracksReport_cube = []
    for Track in Tracks_Filtered:
        views = []
        for view in Track:
            feature = (view[0], view[2][1])
            views.append(feature)
            originTrack_id = view[2][2]
        edges = []
        for Match in Tracks_Matches[originTrack_id]:
            if Match[0] in views and Match[1] in views:
                view1_id = views.index(Match[0])
                view2_id = views.index(Match[1])
                edges.append([view1_id, view2_id])
        TracksReport_cube.append([np.asarray(views), np.asarray(edges)])
    return TracksReport_cube
