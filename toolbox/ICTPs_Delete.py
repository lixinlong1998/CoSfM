import Metashape
import time
import sys

sys.path.append(r'D:\Research\20221223_CoSfM\Release\CFTM_v1.1')
from src_CFTM import ConnectData_Metashape
from src_CFTM import Func_CommonTiePoints
from src_Metashape import FuncMs_CommonTiePoints as MsCTPs

'''Please run this script using the following command in cmd：

metashape.exe -r D:\Research\20221223_CoSfM\Release\CFTM_v1.1\toolbox\ICTPs_Delete.py

'''
#################################################       SETUP      #####################################################
project_path = r"E:\Projects\20230418_CFTM\20240624_Tutorial\Example\cftm_example_project.psx"
chunk_name = ''
epoch_mode = "DATE"  # or "FOLDER"
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
    Points = ConnectData_Metashape.getComponent_Points(chunk, WriteCovariance=False, WriteQuality=False)
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
    print('[Script][TimeCost]    [4]  data analysed:', time.perf_counter() - starttime)

    # [5]  delete MCTPs
    print('[Script]    number of points before delete:', len(chunk.point_cloud.points))
    # ICTPsTypes=['BothRetained', 'Retained1','Retained2','Removed']
    print('[Script]    number of ICTPs:', len(ICTPs_IdList))
    MsCTPs.deleteICTPs(ICTPs_IdList, chunk, Points, TracksEpoch, [])

    # [6] save project cannot update the deleting operation immediately,but it will be updated when reopen the document.
    doc.save()
    print('[time cost]:', time.perf_counter() - starttime0)
