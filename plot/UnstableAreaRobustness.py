# coding:utf-8
import os
import time
import Metashape
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

sys.path.append(r'D:\Research\20221223_CoSfM\Release\CFTM_v1.0')
from src_CFTM import ConnectData_Metashape
from src_Metashape import FuncMs_Basic, FuncMs_Marker
from src_CFTM import Func_ReadData, Func_CommonTiePoints

'''
Please run this script using the following command in cmd:

metashape.exe -r D:/Research/20221223_CoSfM/Release/CFTM_v1.0/plot/UnstableAreaRobustness.py

'''
#################################################       SETUP      #####################################################
# B12 noUR
path_OptmSelect = r'E:\Projects\20230418_CFTM\Baige\Exp6_UnstableMask\Baige_e1e2_Coalign_noCTPs_noUR_ParamMining'
path_Figure1 = r'D:\Research\20221223_CoSfM\Figure\Fig8_unstableMask\20240524_revision/B12_noUR_30_10.svg'
output_path = r'E:\Projects\20230418_CFTM\20240509_Revision\Exp6_UnstableMask'

# B13 noUR
# path_OptmSelect = r'E:\Projects\20230418_CFTM\Baige\Exp6_UnstableMask\Baige_e1e3_Coalign_noCTPs_noUR_ParamMining'
# path_Figure1 = r'D:\Research\20221223_CoSfM\Figure\Fig8_unstableMask\20240524_revision/B13_noUR_30_10.svg'
# output_path = r'E:\Projects\20230418_CFTM\20240509_Revision\Exp6_UnstableMask\B13_NoUA'

iteratNum = 18  # 设置录入的迭代次数 id+1

ax1_ylim = (0, 10)
ax2_ylim = (0, 1)
epoch_mode = 'DATE'
project_path = r"E:\Projects\20230418_CFTM\Baige\Exp5_ParameterMining\Baige_e1e3_Coalign_noCTPs2.psx"
CPsDbs_path = r"E:\Projects\20230418_CFTM\20240509_Revision\Baige_e1e3\B13_CPsdatabase.txt"
bundle_adjustment_args = {
    "fit_f": True,
    "fit_cx": True, "fit_cy": True,
    "fit_b1": True, "fit_b2": True,
    "fit_k1": True, "fit_k2": True,
    "fit_k3": True, "fit_k4": False,
    "fit_p1": True, "fit_p2": True,
    "fit_corrections": False,
    "adaptive_fitting": False,
    "tiepoint_covariance": False
}


#################################################   END OF SETUP   #####################################################
def figure(CTPsQua_List_MixRErr, CPsQua_List_RMSET, iterations, path_Figure1):
    plt.rcParams['font.size'] = 21
    # plot figure
    fig, ax1 = plt.subplots(layout='constrained', figsize=(8, 6), dpi=128)
    ax1.boxplot(CTPsQua_List_MixRErr,
                patch_artist=True,  # 要求用自定义颜色填充盒型图，默认白色填充
                showfliers=False,  # 不显示异常值
                showmeans=True,  # 以点的形式显示均值
                boxprops={'color': 'black', 'facecolor': '#b2b6c5'},
                flierprops={'marker': 'o', 'markerfacecolor': '#e9424a', 'color': 'black'},
                meanprops={'marker': 'o', 'markerfacecolor': '#fcff00', 'markeredgecolor': 'black', 'color': 'black'},
                medianprops={"linestyle": '-', 'color': 'black'})
    ax1.set_ylabel('Epoch Reprojection Error [pix]')

    # ax1.boxplot(CTPsQua_List_RepErr1,
    #             patch_artist=True,  # 要求用自定义颜色填充盒型图，默认白色填充
    #             showfliers=False,  # 不显示异常值
    #             showmeans=True,  # 以点的形式显示均值
    #             boxprops={'color': 'black', 'facecolor': '#a0cac9'},
    #             flierprops={'marker': 'o', 'markerfacecolor': '#e9424a', 'color': 'black'},
    #             meanprops={'marker': 'o', 'markerfacecolor': '#fcff00', 'color': 'black'},
    #             medianprops={"linestyle": '--', 'color': 'black'})
    #
    # ax1.boxplot(CTPsQua_List_RepErr2,
    #             patch_artist=True,  # 要求用自定义颜色填充盒型图，默认白色填充
    #             showfliers=False,  # 不显示异常值
    #             showmeans=True,  # 以点的形式显示均值
    #             boxprops={'color': 'black', 'facecolor': '#a0cac9'},
    #             flierprops={'marker': 'o', 'markerfacecolor': '#e9424a', 'color': 'black'},
    #             meanprops={'marker': 'o', 'markerfacecolor': '#fcff00', 'color': 'black'},
    #             medianprops={"linestyle": '--', 'color': 'black'})

    ax2 = ax1.twinx()
    ax2.plot(iterations, CPsQua_List_RMSET, marker='o', color='red')
    ax2.set_ylabel('Total RMSE of CPs [m]')

    # format
    # 设置横坐标刻度
    plt.xticks(range(1, len(CPsQua_List_RMSET) + 1), iterations)
    # 调整两个y轴的刻度范围
    # ax1.set_xlim(0, iteratNum)
    # ax2.set_xlim(0, iteratNum)
    ax1.set_ylim(ax1_ylim)
    ax2.set_ylim(ax2_ylim)

    # plt.ylim(ylim_Lower, ylim_upper)  # 控制y轴范围在-1.2到1.2之间
    # plt.yticks(list(np.arange(ylim_Lower, ylim_upper, tick_size)) + [ylim_upper])
    # plt.ylabel('Residual[m]', fontsize=14, fontproperties=font)
    # plt.xlabel('Directions', fontsize=14, fontproperties=font)
    # plt.subplots_adjust(left=0.2, bottom=0.2, hspace=0.5)
    # plt.tight_layout()
    # export figure
    plt.savefig(path_Figure1)


def getMixRepErr(CTPsQua_List, iteratNum):
    CTPsQuaIndex_List = []
    for i, CTPsQua in enumerate(CTPsQua_List):
        if i < iteratNum:
            CTPsQuaIndex_List.append([(CTPQua[12] * 2 + CTPQua[13] + CTPQua[14]) / 4 for CTPQua in CTPsQua])
    return CTPsQuaIndex_List


if __name__ == '__main__':
    # 获取CTPs的路径
    CPsRMSET_dict = {}
    CommonTiePoints_path_dict = {}
    for root, dirs, files in os.walk(path_OptmSelect):
        for file in files:
            filename, file_ext = os.path.splitext(file)
            if not file.lower().endswith('.txt'):
                continue
            # 找到CommonTiePoints文件和迭代次数编号
            if "CommonTiePoints" in filename:
                # fileType = filename.split('_')[0]
                iteration_id = int(filename.split('_')[-1])
                CommonTiePoints_path = os.path.join(root, file)
                CommonTiePoints_path_dict[iteration_id] = CommonTiePoints_path

    # 将CTPs导入到project
    doc = Metashape.app.document
    doc.open(project_path)
    chunk = FuncMs_Basic.accessMaxChunk(doc, '')
    camera_ids = ConnectData_Metashape.getCameraIds(chunk)
    cameraPaths = ConnectData_Metashape.getCameraPaths(chunk)
    # assign epochs setting of project
    if epoch_mode == "DATE":
        epochs = Func_CommonTiePoints.getEpoch_byDATE(camera_ids)
    elif epoch_mode == "FOLDER":
        epochs = Func_CommonTiePoints.getEpoch_byFOLDER(cameraPaths)
    CamerasEpoch = Func_CommonTiePoints.analyseCameras_Epoch(camera_ids, cameraPaths, epochs, epoch_mode)
    for iteration_id, CommonTiePoints_path in CommonTiePoints_path_dict.items():
        # 复制chunk
        print('Import CTPs from:', CommonTiePoints_path)
        chunkProcess = chunk.copy(keypoints=False)
        chunkProcess.label = "process"

        # 导入CTPs
        MarkerPoints, MarkerTracks, MarkerInforms = FuncMs_Marker.importMarkersData_Add(CommonTiePoints_path)
        print(f'Add {len(MarkerPoints)} CTPs:')
        FuncMs_Marker.addMarkers(chunkProcess, MarkerPoints, MarkerTracks, 'CTPs')

        # 进行OP
        FuncMs_Basic.bundleAdjustment(chunkProcess, bundle_adjustment_args)

        # 导出CTPs质量
        print('[Script]        Generating CTPs quality report...')
        CTPsData = FuncMs_Marker.getMarkersData_Analyse(chunkProcess, 'CTPs', mode='both')
        CTPsAllQua = FuncMs_Marker.getMarkersAllQua(chunkProcess, CTPsData, camera_ids, CamerasEpoch, 0,
                                                    MarkerFormat='Analyse', Triangulation='linear')
        path_CTPsAllQua = os.path.join(output_path, f'CTPsQua_{iteration_id}.txt')
        path_CTPsAllQuaSta = os.path.join(output_path, f'CTPsQuaSta_{iteration_id}.txt')
        FuncMs_Marker.exportMarkersAllQua(CTPsAllQua, path_CTPsAllQua)
        FuncMs_Marker.reportMarkersAllQua(CTPsAllQua, path_CTPsAllQuaSta)

        # 导出CPs质量
        print('[Script]        Analysing check points...')
        CPsData = FuncMs_Marker.importMarkersData_Analyse(CPsDbs_path, MarkerList=[])
        CPsAllQua = FuncMs_Marker.getMarkersAllQua(chunkProcess, CPsData, camera_ids, CamerasEpoch, 0,
                                                   MarkerFormat='Analyse', Triangulation='nonlinear')
        path_CPsQua = os.path.join(output_path, f'CPsQua_{iteration_id}.txt')
        path_CPsQuaSta = os.path.join(output_path, f'CPsQuaSta_{iteration_id}.txt')
        FuncMs_Marker.exportMarkersAllQua(CPsAllQua, path_CPsQua)
        FuncMs_Marker.reportMarkersAllQua(CPsAllQua, path_CPsQuaSta)

        # 移除临时的chunk
        Metashape.app.document.remove(chunkProcess)
        doc.save()

    # 分析报告
    CTPsQua_List, CTPsQuaSta_List, CPsQuaSta_List = Func_ReadData.readQuas(path_OptmSelect)
    # 读取数据
    CTPsQua_List_ErrPixl = Func_ReadData.getQuaIndList(CTPsQua_List, iteratNum, Ind=12)
    CTPsQua_List_RepErr1 = Func_ReadData.getQuaIndList(CTPsQua_List, iteratNum, Ind=13)
    CTPsQua_List_RepErr2 = Func_ReadData.getQuaIndList(CTPsQua_List, iteratNum, Ind=14)
    CTPsQua_List_MixRErr = getMixRepErr(CTPsQua_List, iteratNum)
    CPsQua_List_RMSET = [CPsQuaSta['Res_T_RMSE(m)'] for i, CPsQuaSta in enumerate(CPsQuaSta_List) if i < iteratNum]
    iterations = range(1, len(CPsQua_List_RMSET) + 1)
    print(CPsQua_List_RMSET)
    print([CPsQua_List_RMSET[i] - CPsQua_List_RMSET[i - 1] for i in range(len(CPsQua_List_RMSET)) if i > 0])
    ####################################################################################################################
    # pre setting
    # plt.rcParams['savefig.facecolor'] = "0.8"
    # plt.rcParams['figure.figsize'] = 4.5, 4.0
    # plt.rcParams['figure.max_open_warning'] = 50
    font = FontProperties()
    font.set_family('Arial')
    # //////////// Fig.1 index boxplot by iterations////////////
    figure(CTPsQua_List_MixRErr, CPsQua_List_RMSET, iterations, path_Figure1)
    # figure(CTPsQua_List_ErrPixl, CPsQua_List_RMSET, iterations, path_Figure2)
    # figure(CTPsQua_List_RepErr1, CPsQua_List_RMSET, iterations, path_Figure3)
    # figure(CTPsQua_List_RepErr2, CPsQua_List_RMSET, iterations, path_Figure4)
