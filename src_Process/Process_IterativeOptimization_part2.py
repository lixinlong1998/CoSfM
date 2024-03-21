from config import *
import Metashape
import random
import csv
import numpy as np
import matplotlib.pyplot as plt
from src_CFTM import Func_Files as FcFile
from src_CFTM import Func_Statistic as FcSta
from src_CFTM import Func_Image as FcImg
from src_Metashape import FuncMs_Marker as MsMarker
from src_Metashape import FuncMs_Basic as FcBasic


def iterater(chunk, args):
    starttime0 = time.perf_counter()
    print('[Script]    Iterative Optimization part 2: run iterations...')

    # [0] unpack arguments
    num_nominated = args["CFTM_args"]["num_nominated"]
    num_selected = args["CFTM_args"]["num_selected"]
    num_max_iterations = args["CFTM_args"]["num_max_iterations"]
    num_inertia = args["CFTM_args"]["num_inertia"]
    threshold_Termination = args["CFTM_args"]["threshold_Termination"]
    unstable_areas_mask_path = args["unstable_areas_mask_path"]
    process_analysis = args["process_args"]["process_analysis"]
    check_points_path = args["check_points_path"]
    iterations_path = args["iterations_path"]
    data_package_path = args["data_package_path"]
    begin_from_breakpoint = args["begin_from_breakpoint"]
    bundle_adjustment_args = args["bundle_adjustment_args"]

    # [1]  导入Markers Grid
    starttime = time.perf_counter()
    print('[Script]    Preparing data...')
    # 创建存储数据的文件
    os.makedirs(iterations_path, exist_ok=True)
    # Import data from metashape
    print('[Script]        Import data from metashape...')
    MarkersGrid = FcFile.readDictionary_Tuple_List(data_package_path + '/MarkersGrid.txt')
    camera_ids = FcFile.readDictionary_Int_Int(data_package_path + '/camera_ids.txt')
    CamerasEpoch = FcFile.readCamerasEpoch(data_package_path + '/CamerasEpoch.txt')
    print('[Script]        data prepared in {0:0.3f} sec.'.format(time.perf_counter() - starttime))

    # [2] 导入稳定区shp以掩膜候选CTPs集合
    if unstable_areas_mask_path:
        MarkersGrid = FcImg.maskCTPs_Unstable(MarkersGrid, unstable_areas_mask_path)
        # 导出MarkersGrid
        path_MarkersMaskbyUnstable = data_package_path + '/MarkersGrid_MaskbyUnstable.txt'
        FcFile.writeMarkersMaskbyUnstable(path_MarkersMaskbyUnstable, MarkersGrid)

    # [3] 创建记录断点和MERE的文件
    path_iterations_MERE = os.path.join(data_package_path, 'iterations_MERE.txt')
    if not os.path.exists(path_iterations_MERE) or not begin_from_breakpoint:
        FcFile.creatIterationsMERE(path_iterations_MERE)

    # [4] 迭代优化
    # create temporary processing chunk
    print('[Script]    Selecting Markers for co-registration...')
    chunkProcess = chunk.copy(keypoints=False)
    chunkProcess.label = "process"
    res_iterations = []
    countDown = 0
    countIncr = 0
    for iteration_id in range(0, num_max_iterations):
        starttime = time.perf_counter()
        if iteration_id == 0:
            MarkersSelected, MarkersAllQua = initialCTPs(chunkProcess, CamerasEpoch, MarkersGrid, camera_ids,
                                                         num_nominated, num_selected)
        else:
            MarkersSelected, MarkersRes = iterateCTPs(chunkProcess, CamerasEpoch, MarkersGrid, camera_ids,
                                                      num_nominated, num_selected, iteration_id)
            # 选出的CTPs添加到原始的chunk中，重新执行BA
            Metashape.app.document.remove(chunkProcess)
            chunkProcess = chunk.copy(keypoints=False)
            chunkProcess.label = "process"
            print('[Script]        Switch back to the original chunk:', chunkProcess.label)

        # 将抽取的点转换为markers添加到chunk
        MarkerGroupName = 'CTP'
        MarkerPoints = [CTPMarker[0] for CTPMarker in MarkersSelected]
        MarkerTracks = [CTPMarker[1] for CTPMarker in MarkersSelected]
        MarkerInforms = [CTPMarker[2] for CTPMarker in MarkersSelected]
        MsMarker.addMarkers(chunkProcess, MarkerPoints, MarkerTracks, MarkerGroupName)
        print(f'[Script]        Add {len(MarkerPoints)} markers into chunk:{chunkProcess.label}.')

        # 执行光束法平差，进行初始化配准
        FcBasic.bundleAdjustment(chunkProcess, bundle_adjustment_args)
        # chunkProcess.optimizeCameras(fit_f=False, fit_cx=True, fit_cy=True, fit_b1=True, fit_b2=True,
        #                              fit_k1=False, fit_k2=True, fit_k3=True, fit_k4=False, fit_p1=False, fit_p2=False,
        #                              fit_corrections=False, adaptive_fitting=False, tiepoint_covariance=False)

        # 导出选中的点
        print('[Script]        Exporting CTPs...')
        path_MarkerSampled = iterations_path + '/CommonTiePoints_{0}.txt'.format(iteration_id)
        MsMarker.exportMarkersData_Add(path_MarkerSampled, MarkerPoints, MarkerTracks, MarkerInforms)

        # 对选中的点进行质量评估，并导出评估报告
        CTPsAllQua = analyseCTPsQuality(chunkProcess, MarkerGroupName, camera_ids, CamerasEpoch, iterations_path,
                                        iteration_id)
        if process_analysis:
            # 绘制Error(pix), RepErr_v3_e1, RepErr_v3_e2的分布和相关性图
            plotCTPsResiduals(CTPsAllQua, iterations_path, iteration_id)

        if process_analysis and check_points_path:
            # 使用外部检查点监测配准效果，并写入文件
            analyseCPs(chunkProcess, check_points_path, camera_ids, CamerasEpoch, iterations_path, iteration_id)

        # check the iteration state
        countDown, countIncr, res_iterations, MixReprojError = iterationState(iteration_id, res_iterations, countDown,
                                                                              countIncr, CTPsAllQua,
                                                                              threshold_Termination, starttime)
        # 打开文本文件，使用 'a' 模式表示追加内容
        with open(path_iterations_MERE, "a") as file:
            fwriter = csv.writer(file, delimiter='\t', lineterminator='\n')
            fwriter.writerow([iteration_id, MixReprojError])

        # check convergence
        if checkConvergence(countDown, countIncr, num_inertia):
            break

    # change the state of iteration recorded in path_iterations_MERE
    changeStateItOp(path_iterations_MERE)
    print('[Script][TimeCost]    Iterative Optimization part 2: run iterations:', time.perf_counter() - starttime0)


def initialCTPs(chunkProcess, CamerasEpoch, MarkersGrid, camera_ids, num_nominated, num_selected):
    '''
    input:
        MarkersGrid = {(r, c):[Marker,...]}
            ,where Marker = [MarkerPoint,track,Information]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
            ,where track = [view,view,...,view]
            ,where Information = [[point_id,Track_id],MarkerKey] or [[CrossMatch,CrossMatch,...,CrossMatch],MarkerKey]
    '''
    # 基于初始条件选择一部分Markers导入到chunk中，进行初始化配准
    print('[Script]    Initializing...')
    print('[Script]        chunk.label:', chunkProcess.label)
    # 基于初始条件选择用于初始配准的点集
    MarkersAllQua = []
    MarkersSelected = []
    MarkersNumber = 0
    GridsNumber = 0
    for grid_id, Markers in MarkersGrid.items():
        MarkersNumber += len(Markers)
        GridsNumber += 1
        # 选取n个点，如果不够则跳过
        if len(Markers) < num_nominated:
            continue
        MarkersAllQua_grid = MsMarker.getMarkersAllQua(chunkProcess, Markers, camera_ids, CamerasEpoch, grid_id,
                                                       MarkerFormat='Add', Triangulation='linear')
        MarkersAllQua += MarkersAllQua_grid
        # 匹配数：越多越好
        # 时相内交会角:越大越好
        MarkersAllQua_array = np.asarray(MarkersAllQua_grid)
        sorted_by_matches = MarkersAllQua_array[np.argsort(MarkersAllQua_array[:, 20])]  # from min to max
        sorted_by_triangl = MarkersAllQua_array[np.argsort(MarkersAllQua_array[:, 18])]  # from min to max
        # 根据排名计算得分
        matches_score = np.asarray([[MarkersAllQua[0], score] for score, MarkersAllQua in enumerate(sorted_by_matches)])
        triangl_score = np.asarray([[MarkersAllQua[0], score] for score, MarkersAllQua in enumerate(sorted_by_triangl)])
        ranked_by_matches_score = matches_score[np.argsort(matches_score[:, 1])]  # from min to max
        ranked_by_triangl_score = triangl_score[np.argsort(triangl_score[:, 1])]  # from min to max
        # 计算得分平均值
        calculate_score = np.asarray([[i, (matches_score + ranked_by_triangl_score[i, 1]) * 0.5] for i, matches_score in
                                      enumerate(ranked_by_matches_score)])
        # 根据平均值进行排名
        Score_ranked = calculate_score[np.argsort(-calculate_score[:, 1])]  # from max to min
        # 抽取的点添加到外部大列表中
        MarkersSelected += [Markers[int(MsMarker.readMarkerKey(i)[1])] for i in Score_ranked[:num_nominated, 0]]
    print('[Script]        Number of enabled grids:', GridsNumber)
    print('[Script]        Number of imported markers:', MarkersNumber)
    print('[Script]        Number of selected markers:', len(MarkersSelected))
    return MarkersSelected, MarkersAllQua


def iterateCTPs(chunkProcess, CamerasEpoch, MarkersGrid, camera_ids, num_nominated, num_selected, iteration_id):
    # 选择一部分点转换为markers添加到chunk中
    print('[Script]    ')
    print('[Script]    Starting loop {0}...'.format(iteration_id))
    print('[Script]        chunk.label:', chunkProcess.label)
    # 基于前一次配准的chunk，在每个网格中计算每个marker的指标，这里采用的指标是MixReprojError
    MarkersRes = {}
    MarkersResList = []
    for grid_id, Markers in MarkersGrid.items():
        MarkersRes_grid = MsMarker.getMarkersRes(chunkProcess, Markers, camera_ids, CamerasEpoch,
                                                 Index='MixReprojError', MarkerFormat='Add', Triangulation='linear')
        MarkersRes[grid_id] = MarkersRes_grid
        MarkersResList += MarkersRes_grid
    # 将指标在全局从小到大排序，前N个的Markers的MarkerKey被记入候选名单
    N = len(MarkersGrid) * num_nominated
    MarkersResList = np.asarray(MarkersResList, dtype=object)
    MarkersResList_ranked = MarkersResList[np.argsort(MarkersResList[:, 1])]  # error_pixel from max to min
    if len(MarkersResList) < N:
        NominatedList = [int(MarkerKey) for MarkerKey in MarkersResList_ranked[:, 0]]
    else:
        NominatedList = [int(MarkerKey) for MarkerKey in MarkersResList_ranked[:N, 0]]
    print('[Script]        Number of nominated markers:', len(NominatedList))
    # 在每个网格中选取k个【指标】最小且记录于NominatedList中的CTPs
    MarkersSelected = []
    for grid_id, MarkersRes_grid in MarkersRes.items():
        # 先看MarkersRes_grid中有多少个点在NominatedList中
        MarkersNominated_Key = [int(Res[0]) for Res in MarkersRes_grid if int(Res[0]) in NominatedList]
        if len(MarkersNominated_Key) <= num_selected:
            MarkersSelected_Key = MarkersNominated_Key
        else:
            # 如果MarkersRes_grid中Nominated Markers的数量大于k个，那么就排序，选最小的k个
            MarkersRes_grid = np.asarray(MarkersRes_grid, dtype=object)
            MarkersRes_grid_ranked = MarkersRes_grid[np.argsort(MarkersRes_grid[:, 1])]  # ErrorPixel from min to max
            MarkersSelected_Key = []
            MarkersSelected_Num = 0
            for i in MarkersRes_grid_ranked[:, 0]:
                if int(i) in MarkersNominated_Key:
                    MarkersSelected_Key.append(int(i))
                    MarkersSelected_Num += 1
                    if MarkersSelected_Num == num_selected:
                        break
        # 将抽取的点添加到外部大列表中
        MarkersSelected += [Marker for Marker in MarkersGrid[grid_id] if int(Marker[2][1]) in MarkersSelected_Key]
    print('[Script]        Number of selected markers:', len(MarkersSelected))
    return MarkersSelected, MarkersRes


def analyseCTPsQuality(chunkProcess, MarkerGroupName, camera_ids, CamerasEpoch, iterations_path, iteration_id):
    # 对选中的点进行质量评估，并导出评估报告
    print('[Script]        Generating CTPs quality report...')
    CTPsData = MsMarker.getMarkersData_Analyse(chunkProcess, MarkerGroupName, mode='both')
    CTPsAllQua = MsMarker.getMarkersAllQua(chunkProcess, CTPsData, camera_ids, CamerasEpoch, 0,
                                           MarkerFormat='Analyse', Triangulation='linear')
    path_CTPsAllQua = iterations_path + '/CTPsQua_{0}.txt'.format(iteration_id)
    path_CTPsAllQuaSta = iterations_path + '/CTPsQuaSta_{0}.txt'.format(iteration_id)
    MsMarker.exportMarkersAllQua(CTPsAllQua, path_CTPsAllQua)
    MsMarker.reportMarkersAllQua(CTPsAllQua, path_CTPsAllQuaSta)
    return CTPsAllQua


def plotCTPsResiduals(CTPsAllQua, iterations_path, iteration_id):
    # 绘制Error(pix), RepErr_v3_e1, RepErr_v3_e2的分布和相关性图
    print('[Script]        Generating distribution of errors...')
    CTPsRes = [[AllQua[12], AllQua[13], AllQua[14]] for AllQua in CTPsAllQua]
    path_Figure_Selected = iterations_path + '/ResidualDistrib_{0}.jpg'.format(iteration_id)
    plot_hist_together(path_Figure_Selected, CTPsRes, bin_width=1)


def analyseCPs(chunkProcess, CPsDbs_path, camera_ids, CamerasEpoch, iterations_path, iteration_id):
    # 使用外部检查点监测配准效果，并写入文件
    print('[Script]        Analysing check points...')
    CPsData = MsMarker.importMarkersData_Analyse(CPsDbs_path, MarkerList=[])
    CPsAllQua = MsMarker.getMarkersAllQua(chunkProcess, CPsData, camera_ids, CamerasEpoch, 0,
                                          MarkerFormat='Analyse', Triangulation='nonlinear')
    path_CPsAllQua = iterations_path + '/CPsQua_{0}.txt'.format(iteration_id)
    path_CPsAllQuaSta = iterations_path + '/CPsQuaSta_{0}.txt'.format(iteration_id)
    MsMarker.exportMarkersAllQua(CPsAllQua, path_CPsAllQua)
    MsMarker.reportMarkersAllQua(CPsAllQua, path_CPsAllQuaSta)


def iterationState(iteration_id, res_iterations, countDown, countIncr, CTPsAllQua, threshold_Termination, starttime):
    # 如果平均error pixel小于0.3，或者res RMSE P增量小于0.001，则停止迭代
    MixReprojError = getResidualIndex(CTPsAllQua)
    print('[Script]    MixReprojError:', MixReprojError)
    if iteration_id == 0:
        res_iterations.append(MixReprojError)
        print('[Script]    convergence countDown:{0}'.format(countDown))
        print('[Script]    loop {0} processed in {1:0.3f} sec:'.format(iteration_id, time.perf_counter() - starttime))
        return countDown, countIncr, res_iterations, MixReprojError
    else:
        det_MixReprojError = res_iterations[-1] - MixReprojError
        res_iterations.append(MixReprojError)

        # 如果迭代正在收敛或发散，我们会进入倒数阶段，连续收敛或发散num_inertia次则迭代停止
        if abs(det_MixReprojError) <= threshold_Termination:
            countDown += 1
            countIncr = 0
            print('[Script]    Note: MixReprojError is nearly converged.')
            print('[Script]    countIncr:{0}'.format(countIncr))
            print('[Script]    countDown:{0}'.format(countDown))
            print(
                '[Script]    loop {0} processed in {1:0.3f} sec:'.format(iteration_id, time.perf_counter() - starttime))
            return countDown, countIncr, res_iterations, MixReprojError
        elif det_MixReprojError < (-1 * threshold_Termination):
            countDown = 0
            countIncr += 1
            print('[Script]    Note: MixReprojError is increasing!')
            print('[Script]    countIncr:{0}'.format(countIncr))
            print('[Script]    countDown:{0}'.format(countDown))
            print(
                '[Script]    loop {0} processed in {1:0.3f} sec:'.format(iteration_id, time.perf_counter() - starttime))
            return countDown, countIncr, res_iterations, MixReprojError
        else:
            countDown = 0
            countIncr = 0
            print('[Script]    countIncr:{0}'.format(countIncr))
            print('[Script]    countDown:{0}'.format(countDown))
            print(
                '[Script]    loop {0} processed in {1:0.3f} sec:'.format(iteration_id, time.perf_counter() - starttime))
            return countDown, countIncr, res_iterations, MixReprojError


def checkConvergence(countDown, countIncr, num_inertia):
    # 是否连续收敛或发散num_inertia次
    if countIncr == num_inertia:
        print(f'[Script]    Stop: MixReprojError increased for {num_inertia} times!')
        return True
    elif countDown == num_inertia:
        print('[Script]    Stop: MixReprojError converged.')
        return True
    else:
        return False


def changeStateItOp(path):
    # 读取文件内容
    with open(path, "r") as file:
        lines = file.readlines()
    if lines:
        lines[0] = "Iteration state: complete\n"
    # 将修改后的内容写回文件
    with open(path, "w") as file:
        file.writelines(lines)


def randomSampling(lst, k):
    """
    在长度为 n 的列表 lst 中均匀抽取 k 个元素，返回包含这 k 个元素的列表和剩余元素的列表
    """
    n = len(lst)
    indices = set()
    while len(indices) < k:
        indices.add(random.randint(0, n - 1))
    return [lst[i] for i in indices], [lst[i] for i in range(n) if i not in indices]


def plot_hist_together(path, CTPsResiduals, bin_width=0.0):
    CTPsResiduals = np.asarray(CTPsResiduals)
    ErrorPix = CTPsResiduals[:, 0]
    RepErrV3e1 = CTPsResiduals[:, 1]
    RepErrV3e2 = CTPsResiduals[:, 2]
    # 计算每个指标的最小值、最大值和范围,根据给定的柱子宽度计算柱子的边界值
    CTPsRes_MIN = np.min(CTPsResiduals, axis=0)
    CTPsRes_MAX = np.max(CTPsResiduals, axis=0)
    x_bins = np.arange(CTPsRes_MIN[0], CTPsRes_MAX[0] + bin_width, bin_width)
    y_bins = np.arange(CTPsRes_MIN[1], CTPsRes_MAX[1] + bin_width, bin_width)
    z_bins = np.arange(CTPsRes_MIN[2], CTPsRes_MAX[2] + bin_width, bin_width)
    # 创建一个大小为6x4英寸的画布和一个Axes对象
    fig, ax = plt.subplots(figsize=(10, 6), dpi=128)
    # 绘制直方图
    ax.hist(ErrorPix, bins=x_bins, alpha=0.3, color='blue', label='Error(pixel)')
    ax.hist(RepErrV3e1, bins=y_bins, alpha=0.3, color='green', label='RepErrV3e1')
    ax.hist(RepErrV3e2, bins=z_bins, alpha=0.3, color='red', label='RepErrV3e2')
    # 添加图例和坐标轴标签
    ax.legend(loc='best')
    ax.set_xlabel('index')
    ax.set_ylabel('frequency')
    # 调整子图之间的距离和周围留白
    plt.subplots_adjust(wspace=0.3, left=0.1, right=0.95, bottom=0.15)
    # 显示图形
    plt.savefig(path)
    plt.close()


def getResidualIndex(MarkersAllQua):
    Distrib_ErrorPixel = FcSta.listStatistic([AllQua[12] for AllQua in MarkersAllQua])
    Distrib_RepErrV3e1 = FcSta.listStatistic([AllQua[13] for AllQua in MarkersAllQua])
    Distrib_RepErrV3e2 = FcSta.listStatistic([AllQua[14] for AllQua in MarkersAllQua])
    MixReprojError = (Distrib_ErrorPixel[0] * 2 + Distrib_RepErrV3e1[0] + Distrib_RepErrV3e2[0]) / 4
    return MixReprojError
