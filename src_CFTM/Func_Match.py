import os

# os.add_dll_directory(r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8\bin')  # for import cv2
# os.add_dll_directory(r'G:\UAV_SOFTWARE\OpenCV\opencv_4_5_0_cuda_11_1_py38\install\x64\vc16\bin')  # for import cv2
import cv2 as cv
import numpy as np

import src_CFTM.Func_Triangulate as Tri
import src_CFTM.Func_Camera as Cam


def innerMatch(features_cube, chunkKeypoints, chunkDescriptors, Cameras, Sensors, camera_ids, criterias, exportReport):
    '''
    input:
        features_cube = [features_cube_camera_id,features_cube_camera_id,...,features_cube_camera_id]
            ,where features_cube_camera_id = [identifier,keypoint_id,keypoint_id,...,keypoint_id]
        chunkKeypoints
        chunkDescriptors
        Cameras
        Sensors
        camera_ids
        criterias = [threshold_RepError, threshold_Angle, *threshold_ProAcc]
    output:
        matches_cube = [matches_cube_pair_id,matches_cube_pair_id,...,matches_cube_pair_id]
            ,where matches_cube_pair_id = [(camera_id1, camera_id2),match,match,...,match]
            ,where match = [[camera_id, projection_info, keypoints_id],[camera_id, projection_info, keypoints_id], matchQuality]
            ,where projection_info = [u,v]
            ,where matchQuality = [distance, point3D]
    '''
    matches_cube = []
    MatchesReport_cube = []
    for i in range(len(features_cube) - 1):
        for j in range(i + 1, len(features_cube)):

            # Convert feature data structure for camera pair of camera_id1 and camera_id2
            Keypoints1, Descriptors1, camera_id1, Keypoints_ids1 = convertFeaturesCubeData(
                features_cube[i], chunkKeypoints, chunkDescriptors, camera_ids)
            Keypoints2, Descriptors2, camera_id2, Keypoints_ids2 = convertFeaturesCubeData(
                features_cube[j], chunkKeypoints, chunkDescriptors, camera_ids)

            # Use brute force matcher and cross checking
            matches_cube_pair_id = bfMatch(Keypoints1, Descriptors1, Keypoints2, Descriptors2,
                                           camera_id1, camera_id2, Keypoints_ids1, Keypoints_ids2)
            # [Exception Handling]: In case there is no match in pair.
            if len(matches_cube_pair_id) == 1:
                # skip this pair
                continue
            # Validate innner matches by two-rays-triangulation
            matches_cube_pair_id_inlier = geometryVerify(matches_cube_pair_id, Cameras, Sensors, criterias)
            # [Exception Handling]: In case there is no match left after two-rays-triangulation verification
            if len(matches_cube_pair_id_inlier) == 1:
                # skip this pair
                continue
            # output valid inner matches
            matches_cube.append(matches_cube_pair_id_inlier)
            if exportReport:
                MatchesReport_cube.append(getMatchesReport(matches_cube_pair_id, matches_cube_pair_id_inlier))
    return matches_cube, MatchesReport_cube


def outerMatch(Track1, Track2, Cameras, Sensors, chunkDescriptors):
    '''
    input:
        Track = [view,view,...view]
            ,where view = [camera_id, projection_info, index_info]
            ,where projection_info = [u, v]
            ,where index_info = [identifier, keypoints_id, originTrack_id]
        chunkDescriptors = 2darray([1darray.shape(128,).dtype(uint8),...])
    output:
        CommonTrackMatches = [CrossMatch,…]
            ,where CrossMatch = [feature1, feature2, MatchQuality]
            ,where feature = [camera_id, keypoint_id]
            ,where MatchQuality = [distance,angles],distance is the L2 normal distance between descriptors of match.
    '''

    # get descriptors of tracks
    def getDescriptors(Track, chunkDescriptors):
        '''return:
        Descriptors = 2darray([1darray.shape(128,).dtype(uint8),...])
        '''
        Descriptors = []
        for view_id, view in enumerate(Track):
            identifier = view[2][0]
            keypoints_id = view[2][1]
            Descriptors.append(chunkDescriptors[identifier][keypoints_id])
        return np.asarray(Descriptors)

    des1 = getDescriptors(Track1, chunkDescriptors)
    des2 = getDescriptors(Track2, chunkDescriptors)

    '''match strategy1：BruteForce Matcher with L2 normal and cross check strategy'''
    bf = cv.BFMatcher(cv.NORM_L2, crossCheck=True)
    # bf.distanceThreshold = 50 # this is default.
    goodMatches = bf.match(des1, des2)

    '''match strategy2：BruteForce Matcher with K Nearest Neighbor and ratio test strategy'''
    # bf = cv.BFMatcher()
    # # [Exception Handling]: In case there is no features could be matched.
    # try:
    #     matches = bf.knnMatch(des1, des2, k=2)
    # except:
    #     return [], []
    # # [Exception notice]: If no good match could be found, it will return empty list
    # goodMatches = []
    # for m, n in matches:
    #     if m.distance < 0.7 * n.distance:  # 0.7 or 0.8
    #         goodMatches.append(m)

    # record triangulate angle
    CommonTrackMatches = []
    for match in goodMatches:
        local_id1 = match.queryIdx
        local_id2 = match.trainIdx
        view1 = Track1[local_id1]
        view2 = Track2[local_id2]
        feature1 = [view1[0], view1[2][1]]
        feature2 = [view2[0], view2[2][1]]
        angle = Tri.TriangulationAngle(Cameras[view1[0]], Cameras[view2[0]],
                                       Tri.Triangulate([view1, view2], Cameras, Sensors), dtype='deg')
        CommonTrackMatches.append([feature1, feature2, [match.distance, angle]])
    return CommonTrackMatches


def bfMatch(key1, des1, key2, des2, camera_id1, camera_id2, key1_ids, key2_ids):
    '''
    input:
        key = 2darray([1darray.shape(6,).dtype(float32),...])
        des = 2darray([1darray.shape(128,).dtype(uint8),...])
        camera_id : index of camera in pair
        key_ids = list[keypoints_id,keypoints_id,...,keypoints_id]
    output:
        matches_pair_id = [(camera_id1, camera_id2),match,match,...,match]
        ,where match = [[camera_id, projection_info, keypoints_id],[camera_id, projection_info, keypoints_id], [distance]]
        ,where projection_info = [u,v]
    '''
    # BruteForce Matcher with L2 normal and cross check strategy
    bf = cv.BFMatcher(cv.NORM_L2, crossCheck=True)
    # [Exception Handling]: In case there is no features could be matched.
    try:
        matches = bf.match(des1, des2)
    except:
        return [(camera_id1, camera_id2)]
    # convert matches data structure
    matches_pair_id = [(camera_id1, camera_id2)]
    for match in matches:
        local_id1 = match.queryIdx
        local_id2 = match.trainIdx
        distance = match.distance
        projection_info1 = [key1[local_id1, 0], key1[local_id1, 1]]
        projection_info2 = [key2[local_id2, 0], key2[local_id2, 1]]
        keypoint_id1 = key1_ids[local_id1]
        keypoint_id2 = key2_ids[local_id2]
        matches_pair_id.append(
            [[camera_id1, projection_info1, keypoint_id1], [camera_id2, projection_info2, keypoint_id2], [distance]])
    return matches_pair_id


def convertFeaturesCubeData(features_cube_camera_id, chunkKeypoints, chunkDescriptors, camera_ids):
    '''
    extract keypoints and descriptors of features by keypoint_ids from chunkDescriptors.
    input:
        features_cube_camera_id = [identifier,keypoint_id,keypoint_id,...,keypoint_id]
    output:
        Keypoints = 2darray([1darray.shape(6,).dtype(float32),...])
        Descriptors = 2darray([1darray.shape(128,).dtype(uint8),...])
        camera_id
    '''
    identifier = features_cube_camera_id[0]
    camera_id = camera_ids[identifier]
    Keypoints_ids = features_cube_camera_id[1:]
    Keypoints = np.asarray([chunkKeypoints[identifier][KPid] for KPid in Keypoints_ids])
    Descriptors = np.asarray([chunkDescriptors[identifier][KPid] for KPid in Keypoints_ids])
    return np.asarray(Keypoints), np.asarray(Descriptors), camera_id, Keypoints_ids


def geometryVerify(matches_pair_id, Cameras, Sensors, criterias):
    '''
    input:
        matches_pair_id = [(camera_id1, camera_id2),match,match,...,match]
            ,where match = [[camera_id, projection_info, keypoints_id],[camera_id, projection_info, keypoints_id], [distance]]
            ,where projection_info = [u,v]
        Cameras
        Sensors
        criterias = [threshold_RepError, threshold_Angle, *threshold_ProAcc]
    output:
        matches_pair_id_verified = [(camera_id1, camera_id2),match,match,...,match]
            ,where match = [[camera_id, projection_info, keypoints_id],[camera_id, projection_info, keypoints_id], matchQuality]
            ,where projection_info = [u,v]
            ,where matchQuality = [distance, point3D]
    '''
    threshold_RepError = criterias[0]
    threshold_Angle = criterias[1]
    camera_id1 = matches_pair_id[0][0]
    camera_id2 = matches_pair_id[0][1]
    Camera_1 = Cameras[camera_id1]
    Camera_2 = Cameras[camera_id2]
    matches_pair_id_verified = [(camera_id1, camera_id2)]
    for match in matches_pair_id[1:]:
        point3D = Tri.Triangulate(match[:-1], Cameras, Sensors)
        # check1: Triangulate successfully?
        if point3D[3] == 0:
            # print('[Script]        Failed: Singular value in triangulate!')
            continue
        # check2: Triangulation angle
        angle = Tri.TriangulationAngle(Camera_1, Camera_2, point3D, dtype='deg')
        if angle < threshold_Angle:
            # print('[Script]        Failed: the triangulation angle is too small!')
            continue
        # check3: Reprojection error
        RepError1 = Cam.ReprojectError(point3D[:-1], Camera_1, Sensors, [match[0][1][0], match[0][1][1]], 'L1')
        RepError2 = Cam.ReprojectError(point3D[:-1], Camera_2, Sensors, [match[1][1][0], match[1][1][1]], 'L1')
        if RepError1 > threshold_RepError or RepError2 > threshold_RepError:
            # print('[Script]        Failed: the reprojection error is too big!')
            continue
        # Add point coordinate information for faster track filter
        match[2].append(point3D)
        matches_pair_id_verified.append(match)
    return matches_pair_id_verified


def getMatchesReport(matches_cube_pair_id, matches_cube_pair_id_inlier):
    '''
    input:
        matches_cube_pair_id = [(camera_id1, camera_id2),match,match,...,match]
            ,where match = [[camera_id, projection_info, keypoints_id],[camera_id, projection_info, keypoints_id], [distance]]
            ,where projection_info = [u,v]
        matches_cube_pair_id_inlier = [(camera_id1, camera_id2),match,match,...,match]
            ,where match = [[camera_id, projection_info, keypoints_id],[camera_id, projection_info, keypoints_id], matchQuality]
            ,where projection_info = [u,v]
            ,where matchQuality = [distance, point3D]
    output:
        MatchesReport_cube_pair_id = [(camera_id1, camera_id2), Matches, inlierMatches]
            Matches = np.array([[keypoint_id, keypoint_id],…])
            inlierMatches = np.array([[keypoint_id, keypoint_id],…])
    '''
    MatchesReport_cube_pair_id = [matches_cube_pair_id[0]]  # [(camera_id1, camera_id2),]
    MatchesReport_cube_pair_matches = []
    for match in matches_cube_pair_id[1:]:
        keypoint1_id = match[0][2]
        keypoint2_id = match[1][2]
        distance = match[2][0]
        MatchesReport_cube_pair_matches.append([keypoint1_id, keypoint2_id, distance])
    MatchesReport_cube_pair_Inliermatches = []
    for match in matches_cube_pair_id_inlier[1:]:
        keypoint1_id = match[0][2]
        keypoint2_id = match[1][2]
        distance = match[2][0]
        MatchesReport_cube_pair_Inliermatches.append([keypoint1_id, keypoint2_id, distance])
    MatchesReport_cube_pair_id.append(np.asarray(MatchesReport_cube_pair_matches))
    MatchesReport_cube_pair_id.append(np.asarray(MatchesReport_cube_pair_Inliermatches))
    return MatchesReport_cube_pair_id
