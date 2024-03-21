from config import *
from src_CFTM import ConnectData_Metashape
from src_CFTM import Func_SpatialDistribution
from src_CFTM import Func_Files
from src_CFTM import Func_CommonTiePoints
from src_CFTM import Func_CTPsGenerator as CTPsGen
from src_Metashape import FuncMs_Camera


def creatTask(chunk, args):
    # unpack arguments
    PJCS_EPSG = args["reference_args"]["PJCS_EPSG"]
    gird_size = args["CFTM_args"]["gird_size_1"]
    for_all_grids = args["for_all_grids"]
    data_package_path = args["data_package_path"]
    epoch_mode = args["epoch_mode"]

    starttime = time.perf_counter()
    print('[Script]    Preparing data...')
    print('[Script]        Import data from metashape...')
    Cameras = ConnectData_Metashape.getComponent_Cameras(chunk, Covariance=True)
    Sensors = ConnectData_Metashape.getComponent_Sensors(chunk)
    Points = ConnectData_Metashape.getComponent_Points(chunk, WriteCovariance=True, WriteQuality=False)
    Tracks = ConnectData_Metashape.getTracks(chunk)
    CoordinateTransform, CoordinateAttribute = ConnectData_Metashape.getComponent_Coordinates(chunk, PJCS_EPSG)
    camera_ids = ConnectData_Metashape.getCameraIds(chunk)
    cameraPaths = ConnectData_Metashape.getCameraPaths(chunk)
    point_ids = ConnectData_Metashape.getPointIds(chunk)
    print('[Script][TimeCost]        Data prepared:', time.perf_counter() - starttime)

    starttime = time.perf_counter()
    print('[Script]    Analysing data...')
    # get sparse cloud boundary
    Boundary_CLCS, Boundary_PJCS = Func_SpatialDistribution.getBoundary(Points)
    # analysis Common Tie Points
    if epoch_mode == "DATE":
        epochs = Func_CommonTiePoints.getEpoch_byDATE(camera_ids)
    elif epoch_mode == "FOLDER":
        epochs = Func_CommonTiePoints.getEpoch_byFOLDER(cameraPaths)
    CamerasEpoch = Func_CommonTiePoints.analyseCameras_Epoch(camera_ids, cameraPaths, epochs, epoch_mode)
    TracksEpoch = Func_CommonTiePoints.analyseTracks_Epoch(Tracks, epochs, CamerasEpoch)
    ICTPs_IdList, ICTPs_Signal = Func_CommonTiePoints.getICTPsIndex(Points, point_ids, TracksEpoch)
    ICTPsCoord_CLCS, ICTPsCoord_PJCS = Func_CommonTiePoints.getICTPsCoord(Points, ICTPs_IdList)
    print('[Script][TimeCost]        Data analysed:', time.perf_counter() - starttime)

    starttime = time.perf_counter()
    print('[Script]    Creating cude data packages...')
    # 1. Creat grid from separate sparse cloud points
    PointsGrid, ICTPsGrid = CTPsGen.buildPointsGrid(Boundary_PJCS, Points, ICTPsCoord_PJCS, gird_size)
    if for_all_grids:
        Cubes = CTPsGen.buildPointsCubes(Boundary_PJCS, Points, TracksEpoch, PointsGrid, gird_size)
    else:
        # 2. Filter out grid without CTPs
        PointsGrid_NoICTPs = CTPsGen.getGrids_NoICTPs(PointsGrid, ICTPsGrid)
        # 3. Construct cube for cell of PointsGrid_NoICTPs
        Cubes = CTPsGen.buildPointsCubes(Boundary_PJCS, Points, TracksEpoch, PointsGrid_NoICTPs, gird_size)
    # 4. Construct cube data packages
    Camera_IdList_cube = {}
    for grid_id, cube in Cubes.items():
        Camera_IdList = FuncMs_Camera.filterCamerasByPoints(PointsGrid[grid_id], Points, Tracks)
        Camera_IdList_cube[grid_id] = Camera_IdList
    print('[Script][TimeCost]        Data packages created:', time.perf_counter() - starttime)

    starttime = time.perf_counter()
    print('[Script]    Saving task...')
    if not os.path.exists(data_package_path):
        os.mkdir(data_package_path)
    Func_Files.writeDictionary_Tuple_List(data_package_path + '/Cubes.txt', Cubes)
    Func_Files.writeDictionary_Tuple_List(data_package_path + '/Cube_CameraIdList.txt', Camera_IdList_cube)
    Func_Files.writeDictionary_Int_List(data_package_path + '/Cameras.txt', Cameras)
    Func_Files.writeDictionary_Int_Int(data_package_path + '/camera_ids.txt', camera_ids)
    Func_Files.writeDictionary_Str_Int(data_package_path + '/cameraPaths.txt', cameraPaths)
    Func_Files.writeCamerasEpoch(data_package_path + '/CamerasEpoch.txt', CamerasEpoch)
    Func_Files.writeDictionary_Int_List(data_package_path + '/Sensors.txt', Sensors)
    Func_Files.writeCoordinate(data_package_path + '/Coordinate.txt', CoordinateTransform, CoordinateAttribute)
    print('[Script][TimeCost]        Dasks saved:', time.perf_counter() - starttime)
