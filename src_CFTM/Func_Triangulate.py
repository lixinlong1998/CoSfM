import numpy as np
import src_CFTM.Func_Camera as Cam
from scipy.optimize import least_squares


def Triangulate(Track, Cameras, Sensors):
    '''
    Multi-View Triangulation(linear)
    input:
        Track = [view,view,...,view]
            ,where view = [camera_id, projection_info]
            ,where projection_info = [u, v, size, track_id]
        Cameras
        Sensors
    output:
        Point_coord_CLCS = np.ndarray([X,Y,Z,1])
    algorithm:
        The linear function is constructed by multi-view direction rays and the linear least squares problem is solved
        by QR decomposition.
    note:
        the coordinate system used in calculation is 'CLCS' & 'LSECS'.
    '''
    # Because of seperation, Track may not have enough image points, but this case should be check before triangulation
    if len(Track) < 2:
        point3D = np.asarray([0, 0, 0, 0])
    else:
        # Construct coefficient matrix A
        A = []
        for view in Track:
            camera_id = view[0]
            projection_info = view[1]
            # acquire camera pose in chunk local coordinate system(CLCS)
            camera_transform_inv = np.linalg.inv(Cameras[camera_id][0][0])
            camera_Rt = np.mat(camera_transform_inv[:-1])
            # transform image point coordinate(in pixel) to the related direction in camera coordinate system
            point_CAM_normal = Cam.cameraUnprojector(projection_info, Cameras[camera_id], Sensors)  # return: [u,v]
            A.append(np.array(camera_Rt[0, :])[0] - np.array(camera_Rt[2, :])[0] * point_CAM_normal[0])
            A.append(np.array(camera_Rt[1, :])[0] - np.array(camera_Rt[2, :])[0] * point_CAM_normal[1])
        A = np.array(A).reshape(len(Track) * 2, 4)
        u, s, vh = np.linalg.svd(A)
        point3D = vh[-1] / vh[-1][-1]
    return point3D  # <class 'numpy.ndarray'> ([X,Y,Z,1])


def Triangulate_N(Point_coord_init, Track, Cameras, Sensors, WeightFunc=2):
    '''
    Multi-View Triangulation(nonlinear)
    input:
        Point_coord_init = np.ndarray([X,Y,Z,1])
        Track = [view,view,...,view]
            ,where view = [camera_id, projection_info]
            ,where projection_info = [u, v, size, track_id]
        Cameras
        Sensors
        WeightFunc: 0,1,2
            0  not use
            1  min(x)/x
            2  e^(-x)
    output:
        Point_coord_esti = np.ndarray([X,Y,Z,1])
    algorithm:
        Reproject the initial coordinates of the point back to the visible camera, calculate the reprojection error,
        and then optimize the point coordinates by minimizing the sum of the squares of all reprojection errors.
    note:
        the coordinate system used in calculation is 'CLCS' & 'LSECS'
    '''
    # Track should contain at least two views
    if len(Track) < 2:
        return np.asarray([0, 0, 0, 0])

    def fun(Point_coord_init):
        if WeightFunc:
            ProjectionsSize = [view[1][2] for view in Track]
            if WeightFunc == 1:
                # weight function1 = min(x)/x = 1/x
                Weights = [min(ProjectionsSize) / ProjectionSize for ProjectionSize in ProjectionsSize]
            elif WeightFunc == 2:
                # weight function2 = e^(-x)
                Weights = [np.exp(-1 * ProjectionSize) for ProjectionSize in ProjectionsSize]
        RepErrors = []
        for view_id, view in enumerate(Track):
            Camera = Cameras[view[0]]
            projection_coord = view[1][:2]
            if WeightFunc:
                w = np.sqrt(Weights[view_id])
                [u_rep, v_rep] = Cam.cameraProjector(Point_coord_init, Camera, Sensors)
                RepError = np.asarray([u_rep - projection_coord[0], v_rep - projection_coord[1]]) * w
            else:
                [u_rep, v_rep] = Cam.cameraProjector(Point_coord_init, Camera, Sensors)
                RepError = np.asarray([u_rep - projection_coord[0], v_rep - projection_coord[1]])
            RepErrors.append(RepError)
        return np.concatenate(RepErrors)

    res = least_squares(fun, Point_coord_init, loss='linear', f_scale=1.0)
    return np.asarray([res.x[0], res.x[1], res.x[2], 1])


def TriangulationAngle(Camera1, Camera2, point3D, dtype='deg'):
    '''
    Triangulation angle
    parameter:  Camera1, Camera2
                point3D = [X,Y,Z,1]
    return: cos of angle
    algorithm: a dot b = |a|*|b|*cos_theta
    note:
    '''
    camera_center1 = Camera1[0][1]  # numpy.ndarray([X, Y, Z])
    camera_center2 = Camera2[0][1]  # numpy.ndarray([X, Y, Z])
    a = camera_center1 - point3D[:-1]
    b = camera_center2 - point3D[:-1]
    # cos_theta goes from 1 to -1, and theta goes from 0 to 180 degree
    cos_theta = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    if dtype == 'rad':
        return cos_theta
    elif dtype == 'deg':
        return np.arccos(cos_theta) * 180 / np.pi


def PointRepErrList(Point3D, Track, Cameras, Sensors):
    '''
    input:
        Point3D = [X,Y,Z,1]
        Track = [view,view,...,view]
            ,where view = [camera_id, projection_info]
            ,where projection_info = [u, v, size, track_id]
        Cameras
        Sensors
    output:
        RepErrList = [RepErr,RepErr,...,RepErr]
            ,where RepErr is in pixel
    this function reproject the given point to its corresponding cameras to calculate the deviations,
    namely 'reproject error'.
    '''
    RepError_List = []
    for view_id, view in enumerate(Track):
        RepError_List.append(
            Cam.ReprojectError(Point3D[:-1], Cameras[view[0]], Sensors, [view[1][0], view[1][1]], 'L1'))
    return RepError_List


def PointRepErrSigma(Point3D, Track, Cameras, Sensors):
    '''
    input:
        Point3D = [X,Y,Z,1]
        Track = [view,view,...,view]
            ,where view = [camera_id, projection_info]
            ,where projection_info = [u, v, size, track_id]
        Cameras
        Sensors
    output:
        RepErrSigma (in pixel)
    this function reproject the given point to its corresponding cameras to calculate the deviations,
    namely 'reproject error', from image observations. the algorithm use the formula:
                                    sqrt(sum(w*v*v)/(n-1))
    ,which result is the standard deviation of bias
    '''
    RepError_square = []
    for view_id, view in enumerate(Track):
        RepError_square.append(
            Cam.ReprojectError(Point3D[:-1], Cameras[view[0]], Sensors, [view[1][0], view[1][1]], 'L2'))
    if len(RepError_square) > 1:
        return np.sqrt(sum(RepError_square) / (len(RepError_square) - 1))
    else:
        return np.sqrt(sum(RepError_square))


def PointTriAngList(Point3D, Track, Cameras):
    '''
    input:
        Point3D = [X,Y,Z,1]
        Track = [view,view,...,view]
            ,where view = [camera_id, projection_info]
            ,where projection_info = [u, v, size, track_id]
        Cameras
    output:
        TriAngList = [TriAng,TriAng,...,TriAng]
            ,where TriAng is in deg
    '''
    TriAng_List = []
    for i in range(len(Track) - 1):
        for j in range(i + 1, len(Track)):
            TriAng_List.append(TriangulationAngle(Cameras[Track[i][0]], Cameras[Track[j][0]], Point3D, dtype='deg'))
    return TriAng_List
