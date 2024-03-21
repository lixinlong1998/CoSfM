from config import *
from src_CFTM import ConnectData_Metashape
from src_CFTM import Func_SpatialDistribution
from src_CFTM import Func_CommonTiePoints as Func_CommonTiePoints
from src_CFTM import Func_CTPsGenerator as CTPsGen
from src_CFTM import Func_Files as FcFile
from src_CFTM import Func_Geometry as FcGeo
from src_Metashape import FuncMs_CommonTiePoints as MsCTPs
from src_Metashape import FuncMs_Marker as MsMarker
from src_Metashape import FuncMs_Basic as FcBasic


def mergeCTPs(chunk, args):
    starttime0 = time.perf_counter()
    print('[Script]    Iterative Optimization part 1: merge CTPs...')

    # [0] unpack arguments
    PJCS_EPSG = args["reference_args"]["PJCS_EPSG"]
    gird_size = args["CFTM_args"]["gird_size_2"]
    for_all_grids = args["for_all_grids"]
    CTP_path = args["CTP_path"]
    reuse_ICTPs = args["reuse_ICTPs"]
    data_package_path = args["data_package_path"]
    epoch_mode = args["epoch_mode"]
    bundle_adjustment_args = args["bundle_adjustment_args"]

    # [1]  Import data from metashape
    starttime = time.perf_counter()
    print('[Script]    Import data from metashape...')
    Points = ConnectData_Metashape.getComponent_Points(chunk, WriteCovariance=True, WriteQuality=False)
    camera_ids = ConnectData_Metashape.getCameraIds(chunk)
    cameraPaths = ConnectData_Metashape.getCameraPaths(chunk)
    point_ids = ConnectData_Metashape.getPointIds(chunk)
    Tracks = ConnectData_Metashape.getTracks(chunk)
    print('[Script][TimeCost]    [1]  data prepared:', time.perf_counter() - starttime)

    # [2]  analysis data
    starttime = time.perf_counter()
    print('[Script]    Analysing data...')
    # get sparse cloud boundary
    Boundary_CLCS, Boundary_PJCS = Func_SpatialDistribution.getBoundary(Points)
    # get ICTPs
    if epoch_mode == "DATE":
        epochs = Func_CommonTiePoints.getEpoch_byDATE(camera_ids)
    elif epoch_mode == "FOLDER":
        epochs = Func_CommonTiePoints.getEpoch_byFOLDER(cameraPaths)
    CamerasEpoch = Func_CommonTiePoints.analyseCameras_Epoch(camera_ids, cameraPaths, epochs, epoch_mode)
    TracksEpoch = Func_CommonTiePoints.analyseTracks_Epoch(Tracks, epochs, CamerasEpoch)
    ICTPs_IdList, ICTPs_Signal = Func_CommonTiePoints.getICTPsIndex(Points, point_ids, TracksEpoch)
    print('[Script][TimeCost]    [2]  data analysed:', time.perf_counter() - starttime)
    # reuse ICTPs
    if for_all_grids and not reuse_ICTPs:
        ICTPsC1_MarkerPoints = []
        ICTPsC1_MarkerTracks = []
        ICTPsC1_MarkerInforms = []
    else:
        # find out the first class of ICTPs, which type name is 'BothRetained'.
        starttime = time.perf_counter()
        print('[Script]    Finding ICTPsC1...')
        ICTPsTyp = MsCTPs.analyseICTPsType(ICTPs_IdList, Points, TracksEpoch)
        ICTPsC1_IdList = [ICTPs_IdList[i] for i, MCTPTyp in enumerate(ICTPsTyp) if MCTPTyp == 'BothRetained']
        # convert track of ICTPsC1 as marker format
        ICTPsC1_MarkerPoints, ICTPsC1_MarkerTracks, ICTPsC1_MarkerInforms = MsMarker.convertPoint2Marker(
            Points, Tracks, ICTPsC1_IdList)
        # export ICTPsC1 marker
        path_ICTPsC1 = data_package_path + '/ICTPsC1.txt'
        MsMarker.exportMarkersData_Add(path_ICTPsC1, ICTPsC1_MarkerPoints, ICTPsC1_MarkerTracks, ICTPsC1_MarkerInforms)
        # ICTPsC1_MarkerPoints, ICTPsC1_MarkerTracks, ICTPsC1_MarkerInforms = MsMarker.importMarkersData_Add(path_ICTPsC1)
        print('[Script][TimeCost]      ICTPsC1 found.:', time.perf_counter() - starttime)

    # [3] read GCTPs
    starttime = time.perf_counter()
    print('[Script]    Reshaping the data...')
    GCTPs_MarkerPoints = []
    GCTPs_MarkerTracks = []
    GCTPs_MarkerInforms = []
    if CTP_path:
        for root, dirs, files in os.walk(CTP_path):
            fileNameList = files
            break
        for fileName in fileNameList:
            grid_id = fileName.split('_')[-1][:-4]
            fileType = fileName.split('_')[0]
            if fileType == 'CommonTiePoint':
                filePath = CTP_path + '/' + fileName
            else:
                continue
            CommonTracks, CommonTracksMatches = FcFile.readCommonTracks(filePath)
            GCTP_MarkerPoints, GCTP_MarkerTracks, GCTP_MarkerInforms = MsMarker.convertCommonTrack2Marker(
                CommonTracks, CommonTracksMatches, chunk)
            GCTPs_MarkerPoints += GCTP_MarkerPoints
            GCTPs_MarkerTracks += GCTP_MarkerTracks
            GCTPs_MarkerInforms += GCTP_MarkerInforms
    # align ICTPsC1 and GCTPs
    MarkerPoints = GCTPs_MarkerPoints + ICTPsC1_MarkerPoints
    MarkerTracks = GCTPs_MarkerTracks + ICTPsC1_MarkerTracks
    MarkerInforms = GCTPs_MarkerInforms + ICTPsC1_MarkerInforms
    print('[Script][TimeCost]    [3]  Align ICTPsC1 and GCTPs.:', time.perf_counter() - starttime)

    # [4] reshape the data structure by new grids
    starttime = time.perf_counter()
    print('[Script]    Reshaping the data...')
    # build new grids
    MarkersGrid = CTPsGen.buildMarkersGrid(Boundary_PJCS, gird_size, MarkerPoints, MarkerTracks, MarkerInforms)
    GridCoords = FcGeo.creatGrid(Boundary_PJCS, gird_size)
    # export MarkersGrid and Grid coordinate
    FcFile.writeDictionary_Tuple_List(data_package_path + '/MarkersGrid.txt', MarkersGrid)
    FcFile.writeDictionary_Tuple_List(data_package_path + '/GridCoords.txt', GridCoords)
    print('[Script][TimeCost]    [4]  MarkersGrid constructed.:', time.perf_counter() - starttime)

    # [5] delete all ICTPs in the chunk
    starttime = time.perf_counter()
    print('[Script]    Deleting all ICTPs...')
    MsCTPs.deleteICTPs(ICTPs_IdList, chunk, Points, TracksEpoch, [])
    # use optimizeCameras() to update the deleting implementation to point cloud in the present process!
    FcBasic.bundleAdjustment(chunk, bundle_adjustment_args)
    print('[Script][TimeCost]    [5]  All ICTPs deleted:', time.perf_counter() - starttime)
    print('[Script][TimeCost]    Iterative Optimization part 1: merge CTPs:', time.perf_counter() - starttime0)
