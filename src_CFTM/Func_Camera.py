import numpy as np

'''
There are typically four types coordinates of a point when projected to the camera,
    1. word coordinate [X,Y,Z] (CLCS are used in this script)
    2. camera coordinate 
            normal      [X_nol,Y_nol,1]
            distorted   [X_dis,Y_dis,1]
    3. image coordinate [u,v]
'''


def BrownDistortModel(point_CAMCS_normal, calibration_params):
    '''
    Brown Distort Model
    [x_nol, y_nol] -> [x_dis, y_dis]
    '''
    # get normal 3D space camera coordinate(nolX,nolY)
    x, y = point_CAMCS_normal
    # camera model parameters(normal model with f,cx,cy,b1,b2,k1,k2,k3,p1,p2)
    k1 = calibration_params[5]
    k2 = calibration_params[6]
    k3 = calibration_params[7]
    p1 = calibration_params[9]
    p2 = calibration_params[10]
    # input the normal point to distort to image coordinate
    r2 = x * x + y * y
    polynomial1 = 1 + k1 * r2 + k2 * r2 ** 2 + k3 * r2 ** 3
    XYdis = [x * polynomial1 + p1 * (r2 + 2 * x * x) + 2 * p2 * x * y,
             y * polynomial1 + p2 * (r2 + 2 * y * y) + 2 * p1 * x * y]
    return XYdis


def BrownCalibrationModel(point_CAMCS_distorted, calibration_params):
    '''
    Brown Calibration Model
    [x_dis, y_dis] -> [x_nol, y_nol]
    '''
    # get normal 3D space camera coordinate(nolX,nolY)
    x, y = point_CAMCS_distorted
    # camera model parameters(normal model with f,cx,cy,b1,b2,k1,k2,k3,p1,p2)
    k1 = calibration_params[5]
    k2 = calibration_params[6]
    k3 = calibration_params[7]
    p1 = calibration_params[9]
    p2 = calibration_params[10]
    # input the normal point to distort to image coordinate
    r2 = x * x + y * y
    polynomial1 = 1 + k1 * r2 + k2 * r2 ** 2 + k3 * r2 ** 3
    XYnol = [(x - p1 * (r2 + 2 * x * x) - 2 * p2 * x * y) / polynomial1,
             (y - p2 * (r2 + 2 * y * y) - 2 * p1 * x * y) / polynomial1]
    return XYnol


def transPoint_CAM2IMA(point_CAM_distorted, calibration_params):
    '''
    [x_dis, y_dis] -> [u, v]
    '''
    x_dis = point_CAM_distorted[0]
    y_dis = point_CAM_distorted[1]
    f = calibration_params[0]
    cx = calibration_params[1]
    cy = calibration_params[2]
    b1 = calibration_params[3]
    b2 = calibration_params[4]
    height = calibration_params[13]
    width = calibration_params[14]
    u = width / 2 + cx + x_dis * (f + b1) + y_dis * b2
    v = height / 2 + cy + y_dis * f
    return [u, v]


def transPoint_IMA2CAM(point_IMA, calibration_params):
    '''
    [u, v] -> [x_dis, y_dis]
    '''
    u = point_IMA[0]
    v = point_IMA[1]
    f = calibration_params[0]
    cx = calibration_params[1]
    cy = calibration_params[2]
    b1 = calibration_params[3]
    b2 = calibration_params[4]
    height = calibration_params[13]
    width = calibration_params[14]
    y_dis = (v - height / 2 - cy) / f
    x_dis = (u - width / 2 - cx - y_dis * b2) / (f + b1)
    return [x_dis, y_dis]


def transPoint_CLCS2CAMCS(Point_coord_CLCS, Camera):
    '''
    [X, Y, Z, *1] -> [x_nol, y_nol]
    '''
    Point_coord_CLCS = np.asarray([Point_coord_CLCS[0], Point_coord_CLCS[1], Point_coord_CLCS[2], 1])
    cameraTransformMatrix = Camera[0][0]  # 4darray
    cameraTransformMatrix_inv = np.linalg.inv(cameraTransformMatrix)
    point_CAMCS = np.matmul(cameraTransformMatrix_inv, Point_coord_CLCS)
    point_CAMCS_normal = point_CAMCS[:2] / point_CAMCS[2]
    return [point_CAMCS_normal[0], point_CAMCS_normal[1]]


def cameraProjector(Point_coord_CLCS, Camera, Sensors):
    '''
    project a 3D point into image coordinate of given photo
    [X, Y, Z, *1] -> [u, v]
    '''
    # get data
    sensors_id = Camera[3]  # int
    calibration_params = Sensors[sensors_id][0]  # list
    # transform a 3D point to CAMCS
    [x_nol, y_nol] = transPoint_CLCS2CAMCS(Point_coord_CLCS, Camera)
    # Brown Distort Model: transform undistorted coordinate to distorted coordinate
    [x_dis, y_dis] = BrownDistortModel([x_nol, y_nol], calibration_params)
    [u, v] = transPoint_CAM2IMA([x_dis, y_dis], calibration_params)
    return [u, v]


def cameraProjector_array(Camera, Sensors, PointList):
    '''
    project a list of 3D points into image coordinate of given photo
    '''
    """Convert 3-D points to 2-D by projecting onto images.
    array_P3D_OBS,(nobservations,3)
    array_CSP_OBS,(nobservations,7)
    array_CMP_OBS,(nobservations,10)
    """
    CSP_qua = array_CSP_OBS[:, :4]
    CSP_rot = Rotation.from_quat(CSP_qua).as_matrix()
    CSP_loc = array_CSP_OBS[:, 4:]

    CSP_rot_inv = np.linalg.inv(CSP_rot)
    Pw_t = array_P3D_OBS - CSP_loc
    P3D_CSP = np.matmul(CSP_rot_inv, Pw_t[:, :, np.newaxis])
    P3D_CSP_normal = P3D_CSP[:, :2] / P3D_CSP[:, 2, np.newaxis]  # 暂时不用翻转图像，到时候可以用unproject方法验证

    nolX = P3D_CSP_normal[:, 0]
    nolY = P3D_CSP_normal[:, 1]

    f = array_CMP_OBS[:, 0, np.newaxis]
    cx = array_CMP_OBS[:, 1, np.newaxis]
    cy = array_CMP_OBS[:, 2, np.newaxis]
    b1 = array_CMP_OBS[:, 3, np.newaxis]
    b2 = array_CMP_OBS[:, 4, np.newaxis]
    k1 = array_CMP_OBS[:, 5, np.newaxis]
    k2 = array_CMP_OBS[:, 6, np.newaxis]
    k3 = array_CMP_OBS[:, 7, np.newaxis]
    p1 = array_CMP_OBS[:, 8, np.newaxis]
    p2 = array_CMP_OBS[:, 9, np.newaxis]

    r2 = np.sum(P3D_CSP_normal ** 2, axis=1)
    x_dis = nolX * (1 + k1 * r2 + k2 * r2 ** 2 + k3 * r2 ** 3) + (p1 * (r2 + 2 * nolX ** 2) + 2 * p2 * nolX * nolY)
    y_dis = nolY * (1 + k1 * r2 + k2 * r2 ** 2 + k3 * r2 ** 3) + (p2 * (r2 + 2 * nolY ** 2) + 2 * p1 * nolX * nolY)
    u = self.width / 2 + cx + x_dis * (
            f + b1) + y_dis * b2  # width是全局变量; u:<class 'numpy.ndarray'>(nobservations, 1)
    v = self.height / 2 + cy + y_dis * f  # height是全局变量; v:<class 'numpy.ndarray'>(nobservations, 1)
    array_PROJ_OBS = np.hstack([u, v])  # <class 'numpy.ndarray'> (nobservations, 2)
    return array_PROJ_OBS


def cameraUnprojector(Projection, Camera, Sensors):
    '''
    project a 2D point into camera coordinate with calibration
    [u, v] -> [x_nol,y_nol]
    Projection = [u,v]
    '''
    # get data
    sensors_id = Camera[3]  # int
    calibration_params = Sensors[sensors_id][0]  # list
    # transform a 2D point to CAMCS
    [x_dis, y_dis] = transPoint_IMA2CAM([Projection[0], Projection[1]], calibration_params)
    # Brown Calibration Model: transform distorted coordinate to undistorted coordinate
    [x_nol, y_nol] = BrownCalibrationModel([x_dis, y_dis], calibration_params)
    return [x_nol, y_nol]


def ReprojectError(Point_coord_CLCS, Camera, Sensors, projection_coord, norm):
    '''
    norm = 'L1' or 'L2'
    reproject a 3D point into image coordinate to calculate the deviation to its image observation.
    [X, Y, Z, *1] + [u, v] -> RE(in pixel)
    '''
    [u_rep, v_rep] = cameraProjector(Point_coord_CLCS, Camera, Sensors)
    reprojection = np.asarray([u_rep - projection_coord[0], v_rep - projection_coord[1]])
    RepErr = np.linalg.norm(reprojection)
    if norm == 'L1':
        return RepErr
    if norm == 'L2':
        return RepErr * RepErr


def getValidProjections(Point_Coord_CLCS, Camera_IdList, Cameras, Sensors):
    '''
    将三维点投影到候选相机中，查看投影点是否在正确的范围内，范围由像片的宽、高和距边界的距离确定。
    output:
        ValidProjections = [ValidProjection,ValidProjection,...,ValidProjection]
            ,where ValidProjection = [u, v, camera_id]
    '''
    ValidProjections = []
    for camera_id in Camera_IdList:
        Camera = Cameras[camera_id]
        w = Sensors[Camera[3]][0][-2]  # 4000
        h = Sensors[Camera[3]][0][-1]  # 6000
        [u, v] = cameraProjector(Point_Coord_CLCS, Camera, Sensors)
        if 0 < u < w and 0 < v < h:
            ValidProjections.append([u, v, camera_id])
    return ValidProjections
