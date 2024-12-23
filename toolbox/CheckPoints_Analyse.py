import time
import sys
import Metashape

sys.path.append(r'D:\Research\20221223_CoSfM\Release\CFTM_v1.2')
from src_CFTM import Func_CommonTiePoints
from src_CFTM import ConnectData_Metashape
from src_Metashape import FuncMs_Marker as MsMarker

'''Please run this script using the following command in cmd：

metashape.exe -r D:\Research\20221223_CoSfM\Release\CFTM_v1.2\toolbox\CheckPoints_Analyse.py

'''
#################################################       SETUP      #####################################################
# project_path = r"E:\Projects\20230418_CFTM\20240624_Tutorial\Example\cftm_example_project.psx"
# CPsDbs_path = r"E:\Projects\20230418_CFTM\20240624_Tutorial\Dataset\CheckPointsDatabase.txt"
CPs_Enabled = []
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
            chunki_images_num = len(chunk_i.cameras)
            if chunki_images_num >= chunk_images_num:
                chunk_images_num = chunki_images_num
                chunk = chunk_i
            else:
                continue

    # [3]  prepare data
    camera_ids = ConnectData_Metashape.getCameraIds(chunk)
    cameraPaths = ConnectData_Metashape.getCameraPaths(chunk)

    # [4]  get CPs
    # assign epochs setting of project
    if epoch_mode == "DATE":
        epochs = Func_CommonTiePoints.getEpoch_byDATE(camera_ids)
    elif epoch_mode == "FOLDER":
        epochs = Func_CommonTiePoints.getEpoch_byFOLDER(cameraPaths)
    CamerasEpoch = Func_CommonTiePoints.analyseCameras_Epoch(camera_ids, cameraPaths, epochs, epoch_mode)

    # [5]  使用外部检查点监测配准效果，并写入文件
    print('[Script]        Analysing check points...')
    CPsData = MsMarker.importMarkersData_Analyse(CPsDbs_path, MarkerList=CPs_Enabled)
    CPsAllQua = MsMarker.getMarkersAllQua(chunk, CPsData, camera_ids, CamerasEpoch, 0,
                                          MarkerFormat='Analyse', Triangulation='nonlinear')
    path_CPsQua = doc.path[0:(len(doc.path) - 4)] + '_CPsQua.txt'
    path_CPsQuaSta = doc.path[0:(len(doc.path) - 4)] + '_CPsQuaSta.txt'
    MsMarker.exportMarkersAllQua(CPsAllQua, path_CPsQua)
    MsMarker.reportMarkersAllQua(CPsAllQua, path_CPsQuaSta)
    print('[Script][TimeCost]    :', time.perf_counter() - starttime0)

    # [6] 显示统计信息
    Report = MsMarker.readMarkersAllQuaReport(path_CPsQuaSta)
    print(f'Project file: {doc.path[0:(len(doc.path) - 4)]}')
    print('\t  P\tZ\tT')
    print('RMSE\t', Report['Res_P_RMSE(m)'], Report['Res_Z_RMSE(m)'], Report['Res_T_RMSE(m)'])
    print('MAE\t', Report['Res_P_MAE(m)'], Report['Res_Z_MAE(m)'], Report['Res_T_MAE(m)'])
    print('AVG\t', Report['Res_P_AVG(m)'], Report['Res_Z_AVG(m)'], Report['Res_T_AVG(m)'])
    print('STD\t', Report['Res_P_STD(m)'], Report['Res_Z_STD(m)'], Report['Res_T_STD(m)'], '\n')
    print('Distrib_ErrorPixel_AVG:\t', Report['Distrib_ErrorPixel_AVG'])
    print('Distrib_RepErr_v3_e1_AVG:\t', Report['Distrib_RepErr_v3_e1_AVG'])
    print('Distrib_RepErr_v3_e2_AVG:\t', Report['Distrib_RepErr_v3_e2_AVG'])
