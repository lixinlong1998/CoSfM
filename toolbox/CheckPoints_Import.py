import os
import sys
import time
import Metashape

sys.path.append(r'D:\Research\20221223_CoSfM\Release\CFTM_v1.1')
from src_CFTM import ConnectData_Metashape
from src_Metashape import FuncMs_Marker as MsMarker

'''Please run this script using the following command in cmd：

metashape.exe -r D:\Research\20221223_CoSfM\Release\CFTM_v1.1\toolbox\CheckPoints_Import.py

Introduction:
给定CPs按MakrerAdd格式的文件地址，将CTPs导入到打开的工程文件中，可以选择是否执行BA，同时还可以导出Metashape为Markers计算的指标
'''
#################################################       SETUP      #####################################################
# e1e2
# project_path = r"E:\20230418_CFTM\20240509_Revision\Baige_e1e2\Projects\B12_CPs.psx"
# CPsDbs_path = r"E:\20230418_CFTM\20240509_Revision\Baige_e1e2\Baige_e1e2_Coalign_CPsdatabase.txt"

# e1e3
# project_path = r"E:\20230418_CFTM\20240509_Revision\Baige_e1e3\Projects\B13_15VGCPs.psx"
# CPsDbs_path = r"E:\20230418_CFTM\20240509_Revision\Baige_e1e3\Projects\B13_20VGCPs_CPsdatabase.txt"

# xiaomojiu e2e3
project_path = r"E:\20230418_CFTM\20240509_Revision\Xiaomojiu_e2e3\Projects\X23_10VGCPs.psx"
CPsDbs_path = r"E:\20230418_CFTM\20240509_Revision\Xiaomojiu_e2e3\X23_15VGCPs_RMSE0.27.txt"

MarkerGroupName = 'CPs'
chunk_name = ''
#################################################   END OF SETUP   #####################################################
# load all markers in each CTPs file to Metashape project.
if __name__ == '__main__':
    starttime0 = time.perf_counter()
    # open document
    if project_path:
        # run script from cmd
        doc = Metashape.app.document
        doc.open(project_path)
    else:
        # run script from GUI
        doc = Metashape.app.document
    doc.save()  # Document.open(): The document is opened in read-only mode because it is already in use.

    # access chunk
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
    camera_ids = ConnectData_Metashape.getCameraIds(chunk)

    # [4]  get CPs
    CPsData = MsMarker.importMarkersData_Analyse(CPsDbs_path, MarkerList=[])

    # [5]  import CPs as Metashape.Markers
    print('Add MarkerPoints:', len(CPsData))
    print('MarkerGroupName:', MarkerGroupName)
    MsMarker.addMarkers_Analyse(chunk, camera_ids, CPsData, MarkerGroupName)

    # 保存文件
    doc.save()
    print('[Script][TimeCost]    :', time.perf_counter() - starttime0)
