import os
import sys
import time
import Metashape

sys.path.append(r'D:\Research\20221223_CoSfM\Release\CFTM_v1.1')
from src_CFTM import ConnectData_Metashape
from src_CFTM import Func_CommonTiePoints
from src_Metashape import FuncMs_CommonTiePoints as MsCTPs
from src_Metashape import FuncMs_Marker as MsMarker

'''Please run this script using the following command in cmd：

metashape.exe -r D:\Research\20221223_CoSfM\Release\CFTM_v1.1\toolbox\ICTPs_Analyse.py

Introduction:
The Common Tie Points(CTPs) analysis consists of three parts:
    Part1: get CTPs by analysis the epoch of views in each track.
    Part2: export CTPs point quality in metashape, which are RU, RE, PA, IC.
    Part3: calculate CTPs quality by script, which are reprojection error, angle, epoch match number.
'''
#################################################       SETUP      #####################################################
project_path = r"E:\Projects\20230418_CFTM\20240624_Tutorial\Example\cftm_example_project.psx"

chunk_name = ''
epoch_mode = "DATE"  # or "FOLDER"
weights = [0, 0, 0, 0, 0, 0, 2, 4, 8, 16]
#################################################   END OF SETUP   #####################################################
if __name__ == '__main__':
    starttime0 = time.perf_counter()
    # [1]  open document
    if project_path:
        # run script from cmd
        doc = Metashape.app.document
        doc.open(project_path)
    else:
        # run script from GUI
        doc = Metashape.app.document

    # [2]  access chunk
    if chunk_name:
        # choose the chunk with given name
        for chunk_i in Metashape.app.document.chunks:
            if chunk_i.label == chunk_name:
                chunk = chunk_i
    else:
        # choose the chunk with maximum images
        chunk_images_num = 0
        for chunk_i in Metashape.app.document.chunks:
            if chunk_i.enabled:
                chunki_images_num = len(chunk_i.cameras)
                if chunki_images_num >= chunk_images_num:
                    chunk_images_num = chunki_images_num
                    chunk = chunk_i
                else:
                    continue

    # [3]  prepare data
    starttime = time.perf_counter()
    print('[Script]    Preparing data...')
    # Import data from metashape
    Cameras = ConnectData_Metashape.getComponent_Cameras(chunk, Covariance=True)
    Sensors = ConnectData_Metashape.getComponent_Sensors(chunk)
    Points = ConnectData_Metashape.getComponent_Points(chunk, WriteCovariance=False, WriteQuality=False)
    # Projections = DataImporting_Metashape.getComponent_Projections(chunk)
    CoordinateTransform, CoordinateAttribute = ConnectData_Metashape.getComponent_Coordinates(chunk, 32648)
    camera_ids = ConnectData_Metashape.getCameraIds(chunk)
    cameraPaths = ConnectData_Metashape.getCameraPaths(chunk)
    point_ids = ConnectData_Metashape.getPointIds(chunk)
    Tracks = ConnectData_Metashape.getTracks(chunk)
    print('[Script][TimeCost]    [3]  data prepared:', time.perf_counter() - starttime)

    # [4]  get Common Tie Points
    starttime = time.perf_counter()
    print('[Script]    Analysing data...')
    # analysis Common Tie Points
    if epoch_mode == "DATE":
        epochs = Func_CommonTiePoints.getEpoch_byDATE(camera_ids)
    elif epoch_mode == "FOLDER":
        epochs = Func_CommonTiePoints.getEpoch_byFOLDER(cameraPaths)
    CamerasEpoch = Func_CommonTiePoints.analyseCameras_Epoch(camera_ids, cameraPaths, epochs, epoch_mode)
    TracksEpoch = Func_CommonTiePoints.analyseTracks_Epoch(Tracks, epochs, CamerasEpoch)
    ICTPs_IdList, ICTPs_Signal = Func_CommonTiePoints.getICTPsIndex(Points, point_ids, TracksEpoch)
    ICTPsTracks, ICTPsTracks_ids = Func_CommonTiePoints.getICTPsTracks(Tracks, Points, epochs, CamerasEpoch,
                                                                       ICTPs_IdList)
    print('[Script][TimeCost]    [4]  data analysed:', time.perf_counter() - starttime)

    # [5]  analysis CTPs
    if not ICTPs_IdList:
        print('no Common Tie Points!')
    else:
        starttime = time.perf_counter()
        print('Common Tie Points:', len(ICTPs_IdList))
        print('[Script]    Analysing CTPs...')
        # analysis the number of each CTP type
        CTPsNum = MsCTPs.analyseICTPsNum(Points, point_ids, TracksEpoch)
        # analysis the spacial distribution of CTPs
        CTPsSpa = MsCTPs.analyseICTPsSpa(Points, weights, ICTPs_IdList)
        # # analysis the type of CTP
        # CTPsTyp = MsCTPs.analyseICTPsType(ICTPs_IdList, Points, TracksEpoch)
        # analysis the density of CTPs
        # CTPsDen = MsCTPs.analyseICTPsDen(Points, ICTPs_IdList)
        # analysis the quality of CTPs
        CTPsQua, CTPsRetriCoord = MsCTPs.analyseICTPsQua(chunk, ICTPsTracks, ICTPs_IdList, Cameras, Sensors)

        ## get the quality of CTPs, which given by metashape
        # CTPsMsQua = MsCTPs.analyseICTPsQua2(chunk, ICTPs_IdList)
        ## construct CTPs Data
        # ICTPsData = MsCTPs.getICTPs(chunk, Points, ICTPs_IdList, ICTPsTracks, CTPsTyp, CTPsDen, CTPsQua, CTPsRetriCoord)
        print('[Script][TimeCost]    [5]  CTPs analysed:', time.perf_counter() - starttime)

        # [6]  generate report and export results
        starttime = time.perf_counter()
        print('[Script]    Exporting CTPs...')
        path_ICTPsQua = doc.path[0:(len(doc.path) - 4)] + '_ICTPsQua.txt'
        path_ICTPsQuaSta = doc.path[0:(len(doc.path) - 4)] + '_ICTPsQuaSta.txt'
        path_ICTPsInfo = doc.path[0:(len(doc.path) - 4)] + '_ICTPsInfo.txt'
        MsMarker.exportMarkersAllQua(CTPsQua, path_ICTPsQua)
        MsCTPs.reportICTPsQua(CTPsQua, path_ICTPsQuaSta)
        MsCTPs.reportICTPsInfo(CTPsNum, CTPsSpa, weights, path_ICTPsInfo)
        print('[Script][TimeCost]    [6]  CTPs exported:', time.perf_counter() - starttime)
    print('[Script][TimeCost]    :', time.perf_counter() - starttime0)
