import Metashape
import numpy as np
from scipy.optimize import least_squares


def MetashapeType2numpyType(array):
    if isinstance(array, Metashape.Vector):
        numpy_array = np.empty([array.size, ], dtype=float)
        for i in range(array.size):
            numpy_array[i] = array[i]
    elif isinstance(array, Metashape.Matrix):
        numpy_array = np.empty([array.size[0], array.size[1]], dtype=float)
        for i in range(array.size[0]):
            for j in range(array.size[1]):
                numpy_array[i, j] = array[i, j]
    else:
        raise Exception("[Script]    parameter should be a Metashape.Matrix or Metashape.Vector!")
    return numpy_array


def Triangulate(Track, cameras):
    '''
    Multi-View Triangulation
    parameter: Track = [[camera_id, projection_info],...,[camera_id, projection_info]]
    return: 3D point coordinate in 3x1 vector
    algorithm: The linear function is constructed by multi-view direction rays and the linear least squares problem is solved by QR decomposition
    note: the coordinate system used in calculation is 'Local Space Euclidean Coordinate System'
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
            projection_coord = Metashape.Vector([projection_info[0], projection_info[1]])
            camera = cameras[camera_id]
            # acquire camera pose in chunk local coordinate system and transform its dtype to numpy matrix
            camera_Rt = np.mat(MetashapeType2numpyType(camera.transform.inv())[:-1])
            # transform image point coordinate(in pixel) to the related direction in camera coordinate system
            normalized_projection = camera.calibration.unproject(projection_coord)  # return: vector(X/Z,Y/Z,1)
            A.append(np.array(camera_Rt[0, :])[0] - np.array(camera_Rt[2, :])[0] * normalized_projection[0])
            A.append(np.array(camera_Rt[1, :])[0] - np.array(camera_Rt[2, :])[0] * normalized_projection[1])
        A = np.array(A).reshape(len(Track) * 2, 4)
        u, s, vh = np.linalg.svd(A)
        point3D = vh[-1] / vh[-1][-1]
    return point3D  # <class 'numpy.ndarray'> ([X,Y,Z,1])


def Triangulate_N(Point_coord_init, Track, cameras, WeightFunc=2):
    '''
    Multi-View Triangulation(nonlinear)
    input:
        Point_coord_init = Metashape.Vector([X,Y,Z,1])
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
        return [0, 0, 0, 0]
    Point_coord_init = Metashape.Vector([Point_coord_init[0], Point_coord_init[1], Point_coord_init[2]])

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
            camera = cameras[view[0]]
            projection_coord = Metashape.Vector([view[1][0], view[1][1]])
            if WeightFunc:
                w = np.sqrt(Weights[view_id])
                reprojection = camera.project(Point_coord_init)
                RepError = (projection_coord - reprojection) * w
            else:
                reprojection = camera.project(Point_coord_init)
                RepError = (projection_coord - reprojection)
            RepErrors.append(MetashapeType2numpyType(RepError))
        return np.concatenate(RepErrors)

    res = least_squares(fun, Point_coord_init, loss='linear', f_scale=1.0)
    return res.x


def TriangulationAngle(camera1, camera2, point3D, dtype='deg'):
    '''
    Triangulation angle
    parameter:  camera1, camera2
                point3D = [X,Y,Z,1]
    return: cos of angle
    algorithm: a dot b = |a|*|b|*cos_theta
    note:
    '''
    Point3D = Metashape.Vector([point3D[0], point3D[1], point3D[2]])
    camera_center1 = MetashapeType2numpyType(camera1.center)  # numpy.ndarray([X, Y, Z])
    camera_center2 = MetashapeType2numpyType(camera2.center)  # numpy.ndarray([X, Y, Z])
    a = camera_center1 - Point3D
    b = camera_center2 - Point3D
    # cos_theta goes from 1 to -1, and theta goes from 0 to 180 degree
    cos_theta = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    # modify cos_theta if it exceed the boundary numerically
    if 1.0 < cos_theta < 1.0 + 1e-6:
        cos_theta = 1.0
    elif -1.0 - 1e-6 < cos_theta < -1.0:
        cos_theta = -1.0
    if dtype == 'rad':
        return cos_theta
    elif dtype == 'deg':
        return np.arccos(cos_theta) * 180 / np.pi


def ReprojectError(Point_coord_CLCS, camera, projection_coord, norm):
    '''
    input:
        Point_coord_CLCS = [X,Y,Z]
        projection_coord = [U,V]
    norm = 'L1' or 'L2'
    reproject a 3D point into image coordinate to calculate the deviation to its image observation.
    [X, Y, Z, *1] + [u, v] -> RE(in pixel)
    '''
    Point_coord_CLCS = Metashape.Vector([Point_coord_CLCS[0], Point_coord_CLCS[1], Point_coord_CLCS[2]])
    residuals = camera.error(Point_coord_CLCS, [projection_coord[0], projection_coord[1]])
    RepErr = np.linalg.norm(MetashapeType2numpyType(residuals))
    if norm == 'L1':
        return RepErr
    if norm == 'L2':
        return RepErr * RepErr


def PointRepErrList(Point_Coord_CLCS, Track, cameras):
    '''
    input:
        Point_Coord_CLCS = [X,Y,Z,1]
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
        camera = cameras[view[0]]
        projection_coord = view[1][0:2]
        RepError_List.append(ReprojectError(Point_Coord_CLCS, camera, projection_coord, 'L1'))
    return RepError_List


def PointRepErrSigma(Point3D, Track, cameras):
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
        camera = cameras[view[0]]
        projection_coord = view[1][0:2]
        RepError_square.append(ReprojectError(Point3D, camera, projection_coord, 'L2'))
    if len(RepError_square) > 1:
        return np.sqrt(sum(RepError_square) / (len(RepError_square) - 1))
    else:
        return np.sqrt(sum(RepError_square))


def PointTriAngList(Point3D, Track, cameras):
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
            camera1 = cameras[Track[i][0]]
            camera2 = cameras[Track[j][0]]
            TriAng_List.append(TriangulationAngle(camera1, camera2, Point3D, dtype='deg'))
    return TriAng_List


def getPointQuality(point, Track, cameras):
    # Criteria 1: reprojections error
    # the average of reprojections error should not exceed reprojection_error_threshold
    RE_list = []
    for camera_id, projection_info in Track:
        camera = cameras[camera_id]
        projection_coord = Metashape.Vector([projection_info[0], projection_info[1]])
        repro_error = camera.error(point, projection_coord)  # Vector([-0.33389, 0.35795])
        RE_list.append(repro_error.norm())
    # Criteria 2: reconstruction uncertainty(triangulation angle)
    # the maximum triangulation angle should not smaller than triangulation_angle_threshold?
    # note: The cosine of theta decreases monotonically from 0 to 180 degree
    RU_list = []
    for i in range(len(Track) - 1):
        for j in range(i + 1, len(Track)):
            camera1 = cameras[Track[i][0]]
            camera2 = cameras[Track[j][0]]
            cos_theta = TriangulationAngle(camera1, camera2, Point3D, dtype='deg')
            angle = np.arccos(cos_theta) * 180 / np.pi  # numpy.arccos() return angle(in rad) from 0 to pi
            RU_list.append(angle)
    # Criteria 3: projection accuracy
    # eliminate the projection with oversize, which means it was detected in blur condition
    PA_list = []
    for camera_id, projection_info in Track:
        try:
            # check the projection.size, which is the projection_info[2]
            PA_list.append(projection_info[2])
        except IndexError:
            # the check points(Metashape.Marker) have no projection.size, so here is for the case
            PA_list.append(1)
    return RE_list, RU_list, PA_list


def assessPointQuality(RE_list, RU_list, PA_list, criterias=[0.3, 5, 8], Print=False):
    reprojection_error_threshold = criterias[0]
    triangulation_angle_threshold = criterias[1]
    projection_accuracy_threshold = criterias[2]
    RE = np.mean(RE_list)
    RU = np.max(RU_list)
    PA = np.mean(PA_list)

    if RE > reprojection_error_threshold or RU < triangulation_angle_threshold or PA > projection_accuracy_threshold:
        if Print:
            if RE > reprojection_error_threshold:
                print('[Script]        Failed: the reprojection error is too big!')
            elif RU < triangulation_angle_threshold:
                print('[Script]        Failed: the triangulation angle is too small!')
            elif PA > projection_accuracy_threshold:
                print('[Script]        Failed: the image observation is too blur!')
        return False
    else:
        if Print:
            print('[Script]        Triangulate successfully.')
            print('[Script]        Mean projection accuracy:', PA)
            print('[Script]        Maximum triangulation angle:', RU)
            print('[Script]        Mean reprojection error:', RE)
        return True
