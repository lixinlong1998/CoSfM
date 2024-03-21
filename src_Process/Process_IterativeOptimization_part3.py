from config import *
from src_CFTM import Func_Files as FcFile
from src_Metashape import FuncMs_Marker as MsMarker
from src_Metashape import FuncMs_Basic as FcBasic


def optimizer(chunk, args):
    starttime0 = time.perf_counter()
    print('[Script]    Iterative Optimization part 3: implement CTPs...')

    # [0] unpack arguments
    project_name = args["project_name"]
    process_analysis = args["process_args"]["process_analysis"]
    check_points_path = args["check_points_path"]
    iterations_path = args["iterations_path"]
    data_package_path = args["data_package_path"]
    bundle_adjustment_args = args["bundle_adjustment_args"]
    MarkerGroupName = 'CTPs'

    # prepare data
    record_MERE = FcFile.readIterationsMERE(os.path.join(data_package_path, 'iterations_MERE.txt'))
    camera_ids = FcFile.readDictionary_Int_Int(os.path.join(data_package_path, 'camera_ids.txt'))
    CamerasEpoch = FcFile.readCamerasEpoch(os.path.join(data_package_path, 'CamerasEpoch.txt'))

    # 根据CTPs的ERE来选择最佳迭代
    iteration_id = min(record_MERE, key=record_MERE.get)
    path_best_CTPs = iterations_path + '/CommonTiePoints_{0}.txt'.format(iteration_id)

    # 导入CTPs
    MarkerPoints, MarkerTracks, MarkerInforms = MsMarker.importMarkersData_Add(path_best_CTPs)
    print('Add MarkerPoints:', len(MarkerPoints))
    print('MarkerGroupName:', MarkerGroupName)
    MsMarker.addMarkers(chunk, MarkerPoints, MarkerTracks, MarkerGroupName)

    # 执行光束法平差
    FcBasic.bundleAdjustment(chunk, bundle_adjustment_args)
    print('[Script][TimeCost]    Iterative Optimization part 3: implement CTPs:', time.perf_counter() - starttime0)
    return chunk
