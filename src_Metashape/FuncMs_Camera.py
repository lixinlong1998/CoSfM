import Metashape
import math
import time
from src_CFTM import Func_Statistic as MsStat
from src_CFTM import Func_Camera as Cam
import src_Metashape.FuncMs_Transform as MsTrans
import src_Metashape.FuncMs_Geometry as MsGeom

'''
access camera attributes
'''


def getCameraIdentifier(camera):
    '''
    make identifier
    '''
    # '2021:04:04 13:40:13' -> 20210404134013
    identifier = ''
    for s in camera.photo.meta['Exif/DateTimeOriginal']:
        if s.isdigit():
            identifier += s
    # 'DSC00009.JPG' -> 09  (# camera.label is DSC00004 in Metashape 1.8.5)
    # identifier += camera.label[-6:-4]
    identifier += camera.label[-2:]
    return int(identifier)


def getCameraErrors(chunk):
    '''
    相机的pose数据与估计值之差，在投影坐标系下
    output:
        CameraErrors[camera_id] = [CameraError,CameraError,...,CameraError]
            ,where CameraError = [CameraError_loc,CameraError_rot]
            ,where CameraError_loc = [camera_loc_Ena, error_metre,
                                      camera_loc_Ref[0], camera_loc_Ref[1], camera_loc_Ref[2],
                                      camera_loc_Est[0], camera_loc_Est[1], camera_loc_Est[2],
                                      camera_loc_Acc[0], camera_loc_Acc[1], camera_loc_Acc[2],
                                      camera_loc_Err[0], camera_loc_Err[1], camera_loc_Err[2]]
            ,where CameraError_rot = [camera_rot_Ena, error_degree,
                                      camera_rot_Ref[0], camera_rot_Ref[1], camera_rot_Ref[2],
                                      camera_rot_Est[0], camera_rot_Est[1], camera_rot_Est[2],
                                      camera_rot_Acc[0], camera_rot_Acc[1], camera_rot_Acc[2],
                                      camera_rot_Err[0], camera_rot_Err[1], camera_rot_Err[2]]
    '''
    cameras = chunk.cameras
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)

    CameraErrors = [-1] * len(cameras)
    camera_reference_enabled = [1] * len(cameras)

    for camera_id, camera in enumerate(cameras):
        if camera.type == Metashape.Camera.Type.Keyframe:
            print('[Script]    the sensor is keyframe!:', camera.label)
            continue  # skipping Keyframes
        elif not camera.transform:
            print('[Script]    camera had no transform!:', camera.label)
            continue  # skipping no transform
        elif not camera.enabled:
            print('[Script]    camera was not enabled!', camera.label)
            continue  # skipping no transform
        elif not camera.reference.enabled:
            camera_reference_enabled[camera_id] = 0
            print('[Script]    camera reference was not enabled:', camera.label)

        # in GGCS(GDCS or PJCS)
        [X, Y, Z, O, P, K] = transCameraPose_MAT2XYZOPK(camera.transform, crs, M)
        camera_loc_Ena = camera.reference.location_enabled
        camera_rot_Ena = camera.reference.rotation_enabled
        camera_loc_Ref = camera.reference.location  # return Vector([X, Y, Z])
        camera_rot_Ref = camera.reference.rotation  # return Vector([A1,A2,A3])
        camera_loc_Est = Metashape.Vector([X, Y, Z])
        camera_rot_Est = Metashape.Vector([O, P, K])
        camera_loc_Acc = camera.reference.location_accuracy  # return Vector([x,y,z]) or None
        camera_rot_Acc = camera.reference.rotation_accuracy  # return Vector([A1,A2,A3]) or None
        camera_loc_Err = camera_loc_Est - camera_loc_Ref  # Vector([dx,dy,dz])
        camera_rot_Err = camera_rot_Est - camera_rot_Ref  # Vector([do,dp,dk])
        error_metre = camera_loc_Err.norm()
        error_degree = camera_rot_Err.norm()
        CameraErrors[camera_id] = [[camera_loc_Ena, error_metre,
                                    camera_loc_Ref[0], camera_loc_Ref[1], camera_loc_Ref[2],
                                    camera_loc_Est[0], camera_loc_Est[1], camera_loc_Est[2],
                                    camera_loc_Acc[0], camera_loc_Acc[1], camera_loc_Acc[2],
                                    camera_loc_Err[0], camera_loc_Err[1], camera_loc_Err[2]],
                                   [camera_rot_Ena, error_degree,
                                    camera_rot_Ref[0], camera_rot_Ref[1], camera_rot_Ref[2],
                                    camera_rot_Est[0], camera_rot_Est[1], camera_rot_Est[2],
                                    camera_rot_Acc[0], camera_rot_Acc[1], camera_rot_Acc[2],
                                    camera_rot_Err[0], camera_rot_Err[1], camera_rot_Err[2]]]
    return CameraErrors


def getCameraError_analysis(CameraErrors):
    CameraErrors_X = []
    CameraErrors_Y = []
    CameraErrors_Z = []
    CameraErrors_XY = []
    CameraErrors_Total = []
    for error in CameraErrors:
        # residual could also be -1, which means camera.reference is not enabled.
        if error == -1:
            continue
        CameraErrors_X.append(error[0])
        CameraErrors_Y.append(error[1])
        CameraErrors_Z.append(error[2])
        CameraErrors_XY.append(math.sqrt(error[0] ** 2 + error[1] ** 2))
        CameraErrors_Total.append(error.norm())
    # [X_RMSE, Y_RMSE, Z_RMSE, XY_RMSE, Total_RMSE,X_MAE, Y_MAE, Z_MAE, XY_MAE, Total_MAE]
    reportCameraError = [MsStat.calculateRMSE(CameraErrors_X),
                         MsStat.calculateRMSE(CameraErrors_Y),
                         MsStat.calculateRMSE(CameraErrors_Z),
                         MsStat.calculateRMSE(CameraErrors_XY),
                         MsStat.calculateRMSE(CameraErrors_Total),
                         MsStat.calculateMAE(CameraErrors_X),
                         MsStat.calculateMAE(CameraErrors_Y),
                         MsStat.calculateMAE(CameraErrors_Z),
                         MsStat.calculateMAE(CameraErrors_XY),
                         MsStat.calculateMAE(CameraErrors_Total)]
    return reportCameraError


def getCameraIdList_Cube(CubesLayers, CubesCorners, chunk, marg):
    '''

    marg = the distance to boundary of image, in pixel
    '''
    # # 必须使用立方体
    # # 这个算法开销太大，对于每一个三维顶点都要遍历所有的相机，暂时放弃
    # cameras = chunk.cameras
    # h = cameras[0].sensor.calibration.height  # VHY
    # w = cameras[0].sensor.calibration.width  # UWX
    # crs = chunk.crs
    # M = chunk.transform.matrix
    #
    # for grid_id, corners in CubesCorners.items():
    #     # corners = [topLeft, topRight, leftBottom, rightBottom, upper, lower]
    #     UpTopLeft = [corners[0][0], corners[0][1], corners[4]]
    #     UpTopRight = [corners[1][0], corners[1][1], corners[4]]
    #     UpLeftBottom = [corners[2][0], corners[2][1], corners[4]]
    #     UpRightBottom = [corners[3][0], corners[3][1], corners[4]]
    #     DownTopLeft = [corners[0][0], corners[0][1], corners[5]]
    #     DownTopRight = [corners[1][0], corners[1][1], corners[5]]
    #     DownLeftBottom = [corners[2][0], corners[2][1], corners[5]]
    #     DownRightBottom = [corners[3][0], corners[3][1], corners[5]]
    #     cube3DVertexes = [UpTopLeft, UpTopRight, UpLeftBottom, UpRightBottom,
    #                       DownTopLeft, DownTopRight, DownLeftBottom, DownRightBottom]
    #     for Vertex3D_PJCS in cube3DVertexes:
    #         Vertex3D_CLCS = MsTrans.transPointCoord_PJCS2CLCS(Vertex3D_PJCS, crs, M)
    #         for camera_id, camera in enumerate(cameras):
    #             # skip invalid camera
    #             if camera.type == Metashape.Camera.Type.Keyframe:
    #                 print('[Script]    the sensor is keyframe!:', camera.label)
    #                 continue  # skipping Keyframes
    #             elif not camera.transform:
    #                 print('[Script]    camera has no transform!:', camera.label)
    #                 continue  # skipping no transform
    #             elif not camera.enabled:
    #                 print('[Script]    camera is not enabled!', camera.label)
    #                 continue  # skipping camera not used
    #             [u, v] = camera.project(Vertex3D_CLCS)
    #             if marg < u < w - marg and marg < v < h - marg:
    #                 CubesLayers[grid_id].append(camera_id)
    return CubesLayers


def getCameraIdList_CubePoints(CubesPoints, Points, Tracks):
    '''
    input:
        CubesPoints = {grid_id:Point_IdList}
    output:
        Camera_IdList_cube = {grid_id:Camera_IdList}
    '''
    Camera_IdList_cube = {}
    for grid_id, Point_IdList in CubesPoints.items():
        if not Point_IdList:
            continue
        Camera_IdList = filterCamerasByPoints(Point_IdList, Points, Tracks)
        Camera_IdList_cube[grid_id] = Camera_IdList
    return Camera_IdList_cube


def transPoint_CLCS2CAMCS(point_coord_CLCS, camera_transform):
    '''
    point_coord_CLCS = [X,Y,Z,1]
    convert the point.coord from CLCS to camera frame
    '''
    point_coord_CAMCSnol = camera_transform.inv() * point_coord_CLCS
    return point_coord_CAMCSnol


def transCameraPose_MAT2XYZOPK(camera_transform, crs, M):
    '''
    convert the camera.transform (in CLCS) to [X,Y,Z,O,P,K](in PJCS & LSECS)
    '''
    camera_center = camera_transform.translation()
    transform = crs.localframe(M.mulp(camera_center)) * M * camera_transform
    xyz = crs.project(M.mulp(camera_center))
    opk = Metashape.Utils.mat2opk(transform.rotation() * Metashape.Matrix.Diag((1, -1, -1)))
    return [xyz[0], xyz[1], xyz[2], opk[0], opk[1], opk[2]]


def transCameraPose_MAT2XYZYPR(camera_transform, crs, M):
    '''
    convert the camera.transform (in CLCS) to [X,Y,Z,Y,P,R](in PJCS & LSECS)
    '''
    camera_center = camera_transform.translation()
    transform = crs.localframe(M.mulp(camera_center)) * M * camera_transform
    xyz = crs.project(M.mulp(camera_center))
    ypr = Metashape.Utils.mat2ypr(transform.rotation() * Metashape.Matrix.Diag((1, -1, -1)))
    return [xyz[0], xyz[1], xyz[2], ypr[0], ypr[1], ypr[2]]


def filterCamerasByPoint(Track):
    '''
    Filter cameras by selected point
    input:
        Track = [view,view,...,view]
    output:
        camera_IdList = [camera_id,camera_id,...,camera_id]
    '''
    Camera_IdList = []
    for view in Track:
        camera_id = view[0]
        Camera_IdList.append(camera_id)
    return Camera_IdList


def filterCamerasByPoints(PointsIdList, Points, Tracks):
    '''
    Filter cameras by selected point
    input:
        pointsIdList = [point_id,point_id,...,point_id]
        point_ids: return point_id by given track_id
        Tracks: return Track by given track_id
    output:
        camera_IdList = [camera_id,camera_id,...,camera_id]
    '''
    Camera_IdList = []
    for point_id in PointsIdList:
        track_id = Points[point_id][6]
        Track = Tracks[track_id]
        Camera_IdList += filterCamerasByPoint(Track)
    return list(set(Camera_IdList))


def filterCamerasByDensePoints(Queries, Camera_IdList_cube, Boundary_PJCS, CamerasEpoch, Cameras, Sensors, crs, M_inv,
                               girdsize, epoch_id):
    '''
    input:
        Queries = {keypoint_id:Query}
            ,where Query = [coordinates,QueryProjections]
            ,where coordinates = np.array(Point_Coord_PJCS) with 3 dim
            ,where QueryProjections = [] it is empty here, but should be [QueryProj,QueryProj,...,QueryProj]
            ,where QueryProj = [u, v, camera_id]
    output:
        ValidQueries = {keypoint_id:Query}
            ,where Query = [coordinates,QueryProjections]
            ,where coordinates = np.array(Point_Coord_PJCS) with 3 dim
            ,where QueryProjections = [QueryProj,QueryProj,...,QueryProj]
            ,where QueryProj = [u, v, camera_id]
    '''
    starttime = time.perf_counter()
    ValidQueries = {}
    for keypoint_id, Query in Queries.items():
        # get value
        Point_Coord_PJCS = Query[0]

        # transPointCoord_PJCS2CLCS
        Point_Coord_PJCS = Metashape.Vector([Point_Coord_PJCS[0], Point_Coord_PJCS[1], Point_Coord_PJCS[2]])
        Point_Coord_CLCS = M_inv.mulp(crs.unproject(Point_Coord_PJCS))

        # 根据坐标所在的grid_id筛选出可视像片
        [r, c] = MsGeom.getGridId(Point_Coord_PJCS, Boundary_PJCS, girdsize)
        Camera_IdList = Camera_IdList_cube[(r, c)]
        Camera_IdList_epoch = [camera_id for camera_id in Camera_IdList if CamerasEpoch[camera_id] == epoch_id]
        if not Camera_IdList_epoch:
            continue

        # 三维坐标投影到候选的相机中，并且判断重投影点是否在像片给定范围内（范围由到图像边界的距离确定）
        ValidProjections = Cam.getValidProjections(Point_Coord_CLCS, Camera_IdList_epoch, Cameras, Sensors)
        if not ValidProjections:
            continue
        else:
            ValidQueries[keypoint_id] = [Point_Coord_PJCS, ValidProjections]
    print('[Script]    [{}]    filter Cameras By DensePoints[Script][TimeCost]    :'.format(epoch_id), time.perf_counter() - starttime)
    return ValidQueries
