import os
import time
import math
import multiprocessing
import sys

sys.path.append(r'G:\AResearchG\20221223_CoSfM\Code')
from src_CFTM import Func_Files
from src_CFTM import ConnectData_COLMAP

'''Please run this script directly.

Introduction:
并行计算每个cubes上的CTPs
--添加断点续跑功能
--输出失败cube
'''
#################################################       SETUP      #####################################################
# //////////////////  set path

# workspace_path = r"I:\20230617_CTPsGenerator_Baige\Baige_e1e2_CFTM/"
# project_name = "Baige_e1e2_CFTM.psx"
# pathList = [('I:/20230418_CTPsGenerator/Adata/FeatureExtraction/e1_20201014/e1_image1_features.db',
#              'I:/20230418_CTPsGenerator/Adata/e1_20201014/image1/'),
#             ('I:/20230418_CTPsGenerator/Adata/FeatureExtraction/e1_20201014/e1_image2_features.db',
#              'I:/20230418_CTPsGenerator/Adata/e1_20201014/image2/'),
#             ('I:/20230418_CTPsGenerator/Adata/FeatureExtraction/e2_20210703/e2_image1_features.db',
#              'I:/20230418_CTPsGenerator/Adata/e2_20210703/image1/'),
#             ('I:/20230418_CTPsGenerator/Adata/FeatureExtraction/e2_20210703/e2_image2_features.db',
#              'I:/20230418_CTPsGenerator/Adata/e2_20210703/image2/')]

workspace_path = "I:/20230617_CTPsGenerator_Baige/Baige_e1e3_CFTM/"
project_name = "Baige_e1e3_CFTM.psx"
pathList = [('I:/20230418_CTPsGenerator/Adata/FeatureExtraction/e1_20201014/e1_image1_features.db',
             'I:/20230418_CTPsGenerator/Adata/e1_20201014/image1/'),
            ('I:/20230418_CTPsGenerator/Adata/FeatureExtraction/e1_20201014/e1_image2_features.db',
             'I:/20230418_CTPsGenerator/Adata/e1_20201014/image2/'),
            ('I:/20230418_CTPsGenerator/Adata/FeatureExtraction/e3_20211203/e3_image1_features.db',
             'I:/20230418_CTPsGenerator/Adata/e3_20211203/image1/'),
            ('I:/20230418_CTPsGenerator/Adata/FeatureExtraction/e3_20211203/e3_image2_features.db',
             'I:/20230418_CTPsGenerator/Adata/e3_20211203/image2/')]

# workspace_path = "L:/20230524_CTPsGenerator_XiaoMoJiu/XiaoMoJiu_e2e3_CFTM/"
# project_name = "XiaoMoJiu_e2e3_CFTM.psx"
# pathList = [('L:/20230524_CTPsGenerator_XiaoMoJiu/XiaoMoJiu_e2e3_CFTM/XiaoMoJiu_e2e3_COLMAP.db',
#              'L:/20230524_CTPsGenerator_XiaoMoJiu/XiaoMoJiu_e2e3_CFTM/images/')]

# default setting
CTP_path = workspace_path + "CommonTiePoints"
DataPackagePath = workspace_path + "ADataPackage"
COLMAPindex = 'image_path'
# COLMAPindex = 'identifier'

# //////////////////  set arguments
'''
girdsize:   The size of the block square, note that subsequent filtering is performed independently for each grid, in meters
poolSize:   The core number of CPU will be used
criterias1: Maximum reprojection error threshold and minimum intersection Angle threshold in feature matching stage, in pixels and degree respectively
criterias2: Maximum reprojection error threshold and minimum intersection Angle threshold in tracks filtering stage, in pixels and degree respectively
ratio_inlier_init:  An initial guess of the average ratio of inliers in track, in the range from 0 to 1
confidence: A probability of being able to get the consensus set by randomly sample certain times, in the range from 0 to 1
threshold_RepError2:    A threshold used for the final assessment of point quality of track, in pixels
threshold_Distance: For track matching, only track pairs whose distance between points is less than a threshold are considered, in meters
'''
girdsize = 64
poolSize = 16
criterias1 = [2, 3]
criterias2 = [1, 3]
ratio_inlier_init = 0.25
confidence = 0.9
threshold_RepError2 = 0.5
threshold_Distance = 3


#################################################   END OF SETUP   #####################################################

def GenerateCTPs(grid_id, corners, Camera_IdList, Cameras, camera_ids, CamerasEpoch, Sensors,
                 CoordinateTransform, CoordinateAttribute, chunkKeypoints, chunkDescriptors, path_folder, arguments):
    # get arguments
    criterias1 = arguments[0]
    criterias2 = arguments[1]
    threshold_RepError2 = arguments[2]
    ratio_inlier_init = arguments[3]
    confidence = arguments[4]
    threshold_Distance = arguments[5] / CoordinateTransform[4]
    showScreen = False
    # creat exporting file path
    exportReport = True  # If enabled, the output file will occupy 7.25 times the original storage space.
    path_CTPs = path_folder + '/CommonTiePoint_{0}.txt'.format(grid_id)
    path_CamerasMask = path_folder + '/CamerasMask_{0}.npy'.format(grid_id)
    path_Matches = path_folder + '/Matches_{0}.npy'.format(grid_id)
    path_TracksReport = path_folder + '/TracksReport_{0}.npy'.format(grid_id)
    path_FilterTracksReport = path_folder + '/FilterTracksReport_{0}.npy'.format(grid_id)
    path_CommonTracksReport = path_folder + '/CommonTracksReport_{0}.npy'.format(grid_id)

    # Step6: Write CommonTracks to file.txt with form of Metashape.Markers
    starttime1 = time.perf_counter()
    Wrif.writeCommonTracks(path_CTPs, CommonTracks, CommonTracksMatches)
    return None


def print_error(value):
    '''
    这个函数可以输出多进程中的报错，但是不会终止多进程
    '''
    print("error: ", value)


# Creat CTPs for each cube by using multi-process
if __name__ == '__main__':
    starttime0 = time.perf_counter()
    arguments = [criterias1, criterias2, threshold_RepError2, ratio_inlier_init, confidence, threshold_Distance]

    # load data
    starttime = time.perf_counter()
    Cubes = Func_Files.readDictionary_Tuple_List(DataPackagePath + '/Cubes.txt')
    Camera_IdList_cube = Func_Files.readDictionary_Tuple_List(DataPackagePath + '/Cube_CameraIdList.txt')
    Cameras = Func_Files.readDictionary_Int_List(DataPackagePath + '/Cameras.txt')
    camera_ids = Func_Files.readDictionary_Int_Int(DataPackagePath + '/camera_ids.txt')
    cameraPaths = Func_Files.readDictionary_Str_Int(DataPackagePath + '/cameraPaths.txt')
    CamerasEpoch = Func_Files.readCamerasEpoch(DataPackagePath + '/CamerasEpoch.txt')
    Sensors = Func_Files.readDictionary_Int_List(DataPackagePath + '/Sensors.txt')
    CoordinateTransform, CoordinateAttribute = Func_Files.readCoordinate(DataPackagePath + '/Coordinate.txt')
    if COLMAPindex == 'identifier':
        chunkKeypoints, chunkDescriptors = ConnectData_COLMAP.getFeatures(pathList, cameraPaths)
    elif COLMAPindex == 'image_path':
        chunkKeypoints, chunkDescriptors = ConnectData_COLMAP.loadFeatures(pathList, cameraPaths)
    else:
        raise Exception("[Script]    COLMAPindex should be 'identifier' or 'image_path'!")

    # Creat output folder
    if not os.path.exists(CTP_path):
        os.mkdir(CTP_path)

    # creat multi-process pool with given process number
    pool = multiprocessing.Pool(processes=poolSize)
    # Divide cubes into batches with number respected to poolSize and adaptive size.
    Grid_ids = [grid_id for grid_id, corners in Cubes.items()]
    batchSize = math.ceil(len(Cubes) / poolSize)
    batches = [Grid_ids[i:i + batchSize] for i in range(0, len(Cubes), batchSize)]
    print('batchSize:', batchSize)
    print('batchs:', len(batches))

    # For each batch, a process pool is started
    results = []
    for batch in batches:
        result = pool.starmap_async(GenerateCTPs,
                                    [(grid_id, Cubes[grid_id], Camera_IdList_cube[grid_id], Cameras, camera_ids,
                                      CamerasEpoch, Sensors, CoordinateTransform, CoordinateAttribute, chunkKeypoints,
                                      chunkDescriptors, CTP_path, arguments)
                                     for grid_id in batch], error_callback=print_error)
        results.append(result)

    # Wait for all process pools to finish executing
    print('Waiting for all subprocesses done...')
    for result in results:
        result.wait()

    # Closing the process pool
    pool.close()
    pool.join()
    print('All processes done!')
    print('[Script][TimeCost]    :', time.perf_counter() - starttime0)
