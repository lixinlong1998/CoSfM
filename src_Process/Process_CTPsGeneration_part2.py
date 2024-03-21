import sys
import math
import json
import multiprocessing
import numpy as np


def readArguments(path):
    with open(path, 'r') as json_file:
        Dictionary = json.load(json_file)
    return Dictionary


def print_error(value):
    '''
    This function can output the error report in the multi-process, but will not terminate the multi-process
    Tip: Before you run a multi-process programme, run a single programme first, then run the multi-process!
    Because the multi-process run time will not report errors, only stuck there without results!
    '''
    print("error: ", value)


# Creat CTPs for each cube by using multi-process
if __name__ == '__main__':

    # get arguments
    ## The first parameter is the name of the script, from the second parameter onwards are additional parameters
    arguments_path = sys.argv[1:][0]
    args = readArguments(arguments_path)
    PATH_CODE = args["code_path"]
    sys.path.append(PATH_CODE)
    from src_CFTM import Func_Files
    from src_CFTM import Func_CTPsGenerator as CTPsGen
    from config import *

    ## unpack arguments
    input_args = [args["CFTM_args"]["criterias_1"],
                  args["CFTM_args"]["criterias_2"],
                  args["CFTM_args"]["threshold_RepError"],
                  args["CFTM_args"]["ratio_inlier_init"],
                  args["CFTM_args"]["confidence"],
                  args["CFTM_args"]["threshold_Distance"]]
    pool_size = args["CFTM_args"]["pool_size"]
    data_package_path = args["data_package_path"]
    CTP_path = args["CTP_path"]
    feature_path = args["feature_path"]
    begin_from_breakpoint = args["begin_from_breakpoint"]
    state_path = args["state_path"]
    state_CFTM_path = args["state_CFTM_path"]

    starttime0 = time.perf_counter()

    # load data
    starttime = time.perf_counter()
    Cubes = Func_Files.readDictionary_Tuple_List(data_package_path + '/Cubes.txt')
    Camera_IdList_cube = Func_Files.readDictionary_Tuple_List(data_package_path + '/Cube_CameraIdList.txt')
    Cameras = Func_Files.readDictionary_Int_List(data_package_path + '/Cameras.txt')
    camera_ids = Func_Files.readDictionary_Int_Int(data_package_path + '/camera_ids.txt')
    cameraPaths = Func_Files.readDictionary_Str_Int(data_package_path + '/cameraPaths.txt')
    CamerasEpoch = Func_Files.readCamerasEpoch(data_package_path + '/CamerasEpoch.txt')
    Sensors = Func_Files.readDictionary_Int_List(data_package_path + '/Sensors.txt')
    CoordinateTransform, CoordinateAttribute = Func_Files.readCoordinate(data_package_path + '/Coordinate.txt')
    # with open(os.path.join(data_package_path, 'chunkKeypoints.json'), 'r') as json_chunkKeypoints:
    #     chunkKeypoints = json.load(json_chunkKeypoints)
    # with open(os.path.join(data_package_path, 'chunkDescriptors.json'), 'r') as json_chunkDescriptors:
    #     chunkDescriptors = json.load(json_chunkDescriptors)
    chunkKeypoints_data = np.load(os.path.join(data_package_path, 'chunkKeypoints.npz'))
    chunkDescriptors_data = np.load(os.path.join(data_package_path, 'chunkDescriptors.npz'))
    chunkKeypoints = {int(key): chunkKeypoints_data[key] for key in chunkKeypoints_data.files}
    chunkDescriptors = {int(key): chunkDescriptors_data[key] for key in chunkDescriptors_data.files}
    print('chunkKeypoints:', len(chunkKeypoints))
    print('chunkDescriptors:', len(chunkDescriptors))

    # Creat output folder
    if not os.path.exists(CTP_path):
        os.mkdir(CTP_path)
    # Creat process state folder
    if not os.path.exists(state_CFTM_path):
        os.mkdir(state_CFTM_path)

    # Creat multi-process pool with given process number
    pool = multiprocessing.Pool(processes=pool_size)
    # Divide cubes into batches with number respected to pool_size and adaptive size.
    # run from the breakpoint?
    if begin_from_breakpoint:
        finished_grids = []
        for file_name in os.listdir(state_CFTM_path):
            file_path = os.path.join(state_CFTM_path, file_name)
            if os.path.isfile(file_path) and file_name.endswith('.txt'):
                grid_id = eval(file_name[len('state_'):-len('.txt')])
                print(type(grid_id))
                finished_grids.append(grid_id)
        Grid_ids = [grid_id for grid_id, corners in Cubes.items() if grid_id not in finished_grids]
    else:
        Grid_ids = [grid_id for grid_id, corners in Cubes.items()]
    # creat batches
    batch_size = math.ceil(len(Cubes) / pool_size)
    batches = [Grid_ids[i:i + batch_size] for i in range(0, len(Cubes), batch_size)]
    print('batch_size:', batch_size)
    print('batchs:', len(batches))

    # For each batch, a process pool is started
    results = []
    for batch in batches:
        result = pool.starmap_async(CTPsGen.GenerateCTPs,
                                    [(grid_id, Cubes[grid_id], Camera_IdList_cube[grid_id], Cameras, camera_ids,
                                      CamerasEpoch, Sensors, CoordinateTransform, CoordinateAttribute, chunkKeypoints,
                                      chunkDescriptors, CTP_path, state_CFTM_path, input_args)
                                     for grid_id in batch], error_callback=print_error)
        results.append(result)

    # Wait for all process pools to finish executing
    print('Waiting for all subprocesses finished...')
    for result in results:
        result.wait()
    # Closing the process pool
    pool.close()
    pool.join()
    print('All task have finished!')

    print('[Script][TimeCost]    :', time.perf_counter() - starttime0)
