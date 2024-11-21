import numpy as np
import time
import sys
import open3d as o3d
import src_CFTM.Func_Camera as Cam
import src_CFTM.Func_Transform as Trans

'''
读取密集点云（外部格式*.ply, *.obj）
重塑点云的存储格式(kd-tree, cross-grids)
查询密集点云点位于的区块对应的cameras
降采样点云
按外部掩膜(*.shp, *.tif)提取密集点云
提取点云关键点
提取点云几何描述符
'''


def GeometryDescriptor(Point_Coord_LCCS, points_selected):
    geoDes = []
    return np.asarray(geoDes)


def extractDPCFeatures111(path, Camera_IdList_cube, Boundary_PJCS, CamerasEpoch, Cameras, Sensors, CrdT, CrdA, offset,
                          geoRadiu, edge, girdsize, epoch_id):
    '''
    input:

    output:
        DPCFeatures = {keypoint_id:feature}
            ,where feature = [coordinates,geoDescriptor,imgDescriptor]
            ,where coordinates = [Point_Coord_LCCS,Point_Coord_PJCS,Point_Coord_CLCS]
            ,where geoDescriptor = []
            ,where imgDescriptor = []
        CamerasQueries = {camera_id:QueryProjections}
            ,where QueryProjections = [QueryProj,QueryProj,...,QueryProj]
            ,where QueryProj = [u, v, keypoint_id, camera_id]
    '''
    starttime = time.perf_counter()
    # 读取密集点云
    DPC = o3d.io.read_point_cloud(path)
    # 构建 KD-Tree
    pcd_tree = o3d.geometry.KDTreeFlann(DPC)
    # 提取ISS关键点坐标
    print('[Script]    [{}]    Extracting ISS keypoints...:'.format(epoch_id))
    isskeypoints = o3d.geometry.keypoint.compute_iss_keypoints(DPC)
    # 对于每个点,提取RKNN点集以计算几何特征描述符，同时将点投影到对应像片中以记录图像特征询问
    CamerasQueries = {camera_id: [] for camera_id, Camera in Cameras.items()}
    DPCFeatures = {}
    print('[Script]    [{}]    ISS keypoints number:', len(isskeypoints.points))
    for keypoint_id in range(len(isskeypoints.points)):
        starttime2 = time.perf_counter()
        starttime1 = time.perf_counter()
        # coordinate
        Point_Coord_LCCS = isskeypoints.points[keypoint_id]
        # search points in radius(RNN), this step also could use KNN or RKNN.
        # k is the number of points searched
        [k, idx, _] = pcd_tree.search_radius_vector_3d(Point_Coord_LCCS, geoRadiu)
        points_selected = np.asarray(DPC.points)[idx, :]
        # geometry descriptor
        # geoDes = GeometryDescriptor(Point_Coord_LCCS, points_selected)
        geoDes = []
        print('[Script]    [{}]    search_radius_vector_3d[Script][TimeCost]    :'.format(epoch_id), time.perf_counter() - starttime1)

        # get georeferenced coordinate
        starttime1 = time.perf_counter()
        Point_Coord_PJCS = Point_Coord_LCCS + offset
        Point_Coord_CLCS = Trans.transPointCoord_PJCS2CLCS(Point_Coord_PJCS, CrdT, CrdA)
        print('[Script]    [{}]    transPointCoord_PJCS2CLCS[Script][TimeCost]    :'.format(epoch_id),
              time.perf_counter() - starttime1)

        # 记录特征点
        DPCFeatures[keypoint_id] = [[Point_Coord_LCCS, Point_Coord_PJCS, Point_Coord_CLCS],
                                    geoDes, []]  # coordinates, geometry descriptor,image descriptor

        # 根据坐标所在的grid筛选出可视像片（重投影点距离像片边界不能超过给定距离）
        starttime1 = time.perf_counter()
        Camera_IdList = Cam.getCorrespondCameras(Point_Coord_PJCS, Camera_IdList_cube, Boundary_PJCS, girdsize)
        print('[Script]    [{}]    Get Correspond Cameras[Script][TimeCost]    :'.format(epoch_id), time.perf_counter() - starttime1)

        # 移除非本时相的相机
        starttime1 = time.perf_counter()
        Camera_IdList_epoch = [camera_id for camera_id in Camera_IdList if CamerasEpoch[camera_id] == epoch_id]
        print('[Script]    [{}]    Get Cameras in epoch[Script][TimeCost]    :'.format(epoch_id), time.perf_counter() - starttime1)

        # 三维坐标投影到候选的相机中，并且判断重投影点是否在像片给定范围内（范围由到图像边界的距离确定）
        starttime1 = time.perf_counter()
        QueryProjections = Cam.getValidProjections(Point_Coord_CLCS, Camera_IdList_epoch, Cameras, Sensors, edge)
        print('[Script]    [{}]    get Valid Projections[Script][TimeCost]    :'.format(epoch_id), time.perf_counter() - starttime1)

        # 关键点通过投影点向像片询问其周围的图像特征，像片先将询问记录下来，待所有关键点询问完毕后，再统一提取特征答复这些询问。
        starttime1 = time.perf_counter()
        for u, v, camera_id in QueryProjections:
            CamerasQueries[camera_id].append([u, v, keypoint_id, camera_id])
        print('[Script]    [{}]    append QueryProjection[Script][TimeCost]    :'.format(epoch_id), time.perf_counter() - starttime1)
        print('[Script]    [{}]    {}[Script][TimeCost]    :'.format(epoch_id, keypoint_id), time.perf_counter() - starttime2)
    print('[Script]    [{}]    [Script][TimeCost]    :'.format(epoch_id), time.perf_counter() - starttime)
    return DPCFeatures, CamerasQueries


def extractDPCFeatures(path, offset, radius_FPFH, num_FPFH, epoch_id):
    '''
    input:

    output:
        DPCFeatures = {keypoint_id:DPCFeature}
            ,where feature = [coordinates,geoDescriptor,imgDescriptor]
            ,where coordinates = np.array(Point_Coord_LCCS) with 3 dim
            ,where geoDescriptor = np.array(geoDescriptor) with 33 dim
            ,where imgDescriptor = []
        Queries = {keypoint_id:Query}
            ,where Query = [coordinates,QueryProjections]
            ,where coordinates = np.array(Point_Coord_PJCS) with 3 dim
            ,where QueryProjections = [] it is empty here, but should be [QueryProj,QueryProj,...,QueryProj]
            ,where QueryProj = [u, v, camera_id]
    '''
    starttime = time.perf_counter()
    # 读取密集点云
    DPC = o3d.io.read_point_cloud(path)
    DPC = DPC.voxel_down_sample(voxel_size=0.1)
    # 构建 KD-Tree
    pcd_tree = o3d.geometry.KDTreeFlann(DPC)

    # 提取ISS关键点坐标
    tic = time.perf_counter()
    print('[Script]    [{}]    Extracting ISS keypoints...:'.format(epoch_id))

    isskeypoints = o3d.geometry.keypoint.compute_iss_keypoints(DPC)
    print('[Script]    [{}]    [Script][TimeCost]    :'.format(epoch_id), time.perf_counter() - tic)

    # 计算FPFH描述符，33维向量
    print('[Script]    [{}]    Compute FPFH descriptors...:'.format(epoch_id))
    pcd_fpfh = o3d.pipelines.registration.compute_fpfh_feature(
        DPC, o3d.geometry.KDTreeSearchParamHybrid(radius=radius_FPFH, max_nn=num_FPFH))  # 2m,100

    # DPCFeatures, Queries
    DPCFeatures = {}
    Queries = {}
    for keypoint_id in range(len(isskeypoints.points)):
        # coordinate
        Point_Coord_LCCS = isskeypoints.points[keypoint_id]
        # get georeferenced coordinate
        Point_Coord_PJCS = Point_Coord_LCCS + offset
        # 查询距离最近的点
        [_, idx, _] = pcd_tree.search_knn_vector_3d(isskeypoints.points[keypoint_id], 1)
        # 获取对应点的FPFH描述符
        geoDescriptor = pcd_fpfh.data[:, idx[0]]  # numpy 的一维数组，长度为 33
        DPCFeatures[keypoint_id] = [Point_Coord_LCCS, geoDescriptor, []]
        Queries[keypoint_id] = [Point_Coord_PJCS, []]
    print('[Script]    [{}]    [Script][TimeCost]    :'.format(epoch_id), time.perf_counter() - starttime)
    return DPCFeatures, Queries


def ImageDescriptor(CamerasQueries, ):
    # 在matlab中根据清单提取HOPC特征
    for camera_id, QuiryList in CameraQuiryList.items():
        # 读取 .mat 文件
        mat = scipy.io.loadmat(
            'G:/AResearchG/20221223_CoSfM/Experiments/HOPC-descriptor/images/mydata.mat'.format(camera_id))
        Des = mat['Des']  # 将 '数组变量名' 替换为实际的变量名
        for i in range(Des.shape(0) / 72):
            descriptor = Des[i]
            keypoint_id = QuiryList[i][2]
            KeypointsDescriptor[keypoint_id][0].append(descriptor)
