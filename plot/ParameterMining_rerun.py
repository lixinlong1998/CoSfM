import os
import time
import Metashape
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.cm import get_cmap
from matplotlib.colors import LinearSegmentedColormap

sys.path.append(r'D:\Research\20221223_CoSfM\Release\CFTM_v1.0')
from src_CFTM import ConnectData_Metashape
from src_Metashape import FuncMs_Basic, FuncMs_Marker
from src_CFTM import Func_ReadData, Func_CommonTiePoints

'''20240523更新
更新了CPs后，计算需要重新进行；
根据Baige_e1e2_Coalign_noCTPs_noUR_ParamMining_15CPs中每个参数组合中精度最高的interation对应的CTPs重新导入到project中进行OP
然后基于新的CPs更新RMSE，并将结果输出为图片.

Please run this script using the following command in cmd：

metashape.exe -r D:\Research\20221223_CoSfM\Release\CFTM_v1.0\plot\ParameterMining.py

'''

#################################################       SETUP      #####################################################
path_OptmSelect = r'E:\Projects\20230418_CFTM\Baige\Exp5_ParameterMining\Baige_e1e2_Coalign_noCTPs_noUR_ParamMining_15CPs/'
path_Figure1 = r'E:\Projects\20230418_CFTM\Baige\Exp5_ParameterMining\Baige_e1e2_Coalign_noCTPs_noUR_ParamMining\ParamMining.png'
List_lowerBound = [10, 20, 30, 40, 50, 60]
List_selectNumber = [1, 5, 10, 15, 20, 25]

project_path = r"E:\Projects\20230418_CFTM\Baige\Exp5_ParameterMining\Baige_e1e2_Coalign_noCTPs2.psx"
CPsDbs_path = r"E:\Projects\20230418_CFTM\20240509_Revision\Baige_e1e2\B12_CPsdatabase.txt"
output_path = r'E:\Projects\20230418_CFTM\20240509_Revision\Exp5_ParameterMining'

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
epoch_mode = 'DATE'


#################################################   END OF SETUP   #####################################################

def figure_surface(X, Y, Z, path):
    # 绘制三维柱状图
    fig = plt.figure()

    # dx = dy = 0.5 * np.ones_like(zpos)
    # dz = hist.ravel()

    # ax = fig.add_subplot(111, projection='3d')
    # ax.bar3d(X.ravel(), Y.ravel(), np.zeros_like(Z).ravel(), 5, 3, Z.ravel())
    # ax.set_zlim(0, 0.1)

    ax = fig.add_subplot(111, projection='3d')
    # ax.plot_wireframe(X, Y, Z,color='black')
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='black', vmax=0.4, vmin=0.19)  # cmap='viridis'

    # 添加颜色条
    cbar = plt.colorbar(surf, shrink=0.6, aspect=10)
    # ticks = [0.15, 0.16, 0.17, 0.18, 0.19]
    # ticklabel = [str(i) for i in ticks]
    # cbar.set_ticks(ticks)
    # cbar.ax.set_yticklabels(ticklabel, fontsize=14)

    # 设置坐标轴标签
    ax.set_xticks(x)
    ax.set_yticks(y)
    # ax.set_zticks([0, 0.05, 0.1])
    ax.set_xlabel('Nominated number')
    ax.set_ylabel('Selected number')
    ax.set_zlabel('RMSE of CPs [m]')
    # 显示图形
    plt.savefig(path)
    plt.show()


def figure_bar(X, Y, Z, path):
    plt.rcParams['font.size'] = 14  # 设置全局字体大小
    # 绘制三维柱状图
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.bar3d(X.ravel() - 4, Y.ravel() - 2, np.zeros_like(Z).ravel(), 8, 4, Z.ravel(), edgecolor='black', color='white',
             shade=True)
    # ax.set_zlim(0.19, 0.5)

    # dx = dy = 0.5 * np.ones_like(zpos)
    # dz = hist.ravel()
    # # 创建颜色映射
    # cmap = get_cmap('viridis')
    # for r in range(Z.shape[0]):
    #     for c in range(Z.shape[1]):
    #         x = X[r, c] + 0.5  # 柱子中心位于格点
    #         y = Y[r, c] + 0.5  # 柱子中心位于格点
    #         z = Z[r, c]
    #         color = cmap(z)  # 根据高度值获取颜色
    #         ax.bar3d(x, y, 0, 4, 4, z, color=color, edgecolor='black')

    # 设置坐标轴标签
    ax.set_xticks(x)
    ax.set_yticks(y)
    # ax.set_zticks([0, 0.05, 0.1])
    ax.set_xlabel('Nominated number')
    ax.set_ylabel('Selected number')
    ax.set_zlabel('Total RMSE of CPs [m]')
    # 显示图形
    plt.savefig(path)
    plt.show()


if __name__ == '__main__':
    CPsRMSET_ParaMining = {}
    CommonTiePoints_path_dirt = {}
    for root, dirs, files in os.walk(path_OptmSelect):
        folderNameList = dirs
        for folderName in folderNameList:
            # 整合test_id文件中的CPsQuaSta和CTPsQuaSta
            test_id = eval(folderName.split('_')[-1])
            print(type(test_id))
            print('test_id:', test_id)
            folderPath = root + '/' + folderName
            # CPsCTPsQuaSta_db: rows:interation, cols:merged QuaSta
            # CPsCTPsQuaSta_db_ids: the iteration but with order respect to system files, such as 1,10,11,12,2,21,3,4..
            CPsCTPsQuaSta_db, CPsCTPsQuaSta_db_ids = Func_ReadData.mergeCPsCTPsQuaSta(root, [folderName])

            # 提取对应的数据列
            rank_col_id = 8
            db = np.asarray(CPsCTPsQuaSta_db)
            db_ids = np.asarray(CPsCTPsQuaSta_db_ids)
            CPsRMSET = list(db[np.argsort(db_ids)][:, rank_col_id])
            # print(CPsRMSET)

            # 找到最小RMSET对应的iteration得到的CTPs
            Min_CPsRMSET_Interation_id = CPsRMSET.index(min(CPsRMSET))
            CommonTiePoints_path = os.path.join(root, folderName, f'CommonTiePoints_{Min_CPsRMSET_Interation_id}.txt')
            CommonTiePoints_path_dirt[test_id] = CommonTiePoints_path
            # print('min(CPsRMSET):', min(CPsRMSET))
            # print('CommonTiePoints_path:', CommonTiePoints_path)

            # 显示ERROR_Pixel
            # db_ranked = db[np.argsort(db[:, rank_col_id])]  # Res_T_RMSE from min to max
            # db_ids_ranked = db_ids[np.argsort(db[:, rank_col_id])]
            # ERROR_Pixel0 = db[np.argsort(db_ids)][:, 101]
            # ERROR_Pixel1 = db[np.argsort(db_ids)][:, 113]
            # ERROR_Pixel2 = db[np.argsort(db_ids)][:, 125]
            # ERROR_Pixel = (ERROR_Pixel0 * 2 + ERROR_Pixel1 + ERROR_Pixel2) / 4
            # print(ERROR_Pixel)
        break

    # 将CTPs导入到project中
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
    for test_id, CommonTiePoints_path in CommonTiePoints_path_dirt.items():
        # 创建临时chunk
        print('Import CTPs from:', CommonTiePoints_path)
        chunkProcess = chunk.copy(keypoints=False)
        chunkProcess.label = "process"

        # 导入CTPs
        MarkerPoints, MarkerTracks, MarkerInforms = FuncMs_Marker.importMarkersData_Add(CommonTiePoints_path)
        print(f'Add {len(MarkerPoints)} CTPs:')
        FuncMs_Marker.addMarkers(chunkProcess, MarkerPoints, MarkerTracks, 'CTPs')

        # 运行OP
        FuncMs_Basic.bundleAdjustment(chunkProcess, bundle_adjustment_args)

        # 使用新的CPs评估RMSE
        print('[Script]        Analysing check points...')
        CPsData = FuncMs_Marker.importMarkersData_Analyse(CPsDbs_path, MarkerList=[])
        CPsAllQua = FuncMs_Marker.getMarkersAllQua(chunkProcess, CPsData, camera_ids, CamerasEpoch, 0,
                                                   MarkerFormat='Analyse', Triangulation='nonlinear')
        path_CPsQua = os.path.join(output_path, f'{test_id}_CPsQua.txt')
        path_CPsQuaSta = os.path.join(output_path, f'{test_id}_CPsQuaSta.txt')
        FuncMs_Marker.exportMarkersAllQua(CPsAllQua, path_CPsQua)
        FuncMs_Marker.reportMarkersAllQua(CPsAllQua, path_CPsQuaSta)

        # 记录RMSE
        Report = FuncMs_Marker.readMarkersAllQuaReport(path_CPsQuaSta)
        CPsRMSET_ParaMining[test_id] = Report['Res_T_RMSE(m)']

        Metashape.app.document.remove(chunkProcess)
        doc.save()

    print(CPsRMSET_ParaMining)
    CPsRMSET_PM_ListK = [CPsRMSET for test_id, CPsRMSET in CPsRMSET_ParaMining.items()]
    CPsRMSET_PM_ListV = [test_id for test_id, CPsRMSET in CPsRMSET_ParaMining.items()]
    print(max(CPsRMSET_PM_ListK))
    print(CPsRMSET_PM_ListV[CPsRMSET_PM_ListK.index(max(CPsRMSET_PM_ListK))])
    print(min(CPsRMSET_PM_ListK))
    print(CPsRMSET_PM_ListV[CPsRMSET_PM_ListK.index(min(CPsRMSET_PM_ListK))])

    # 定义数据
    x = np.asarray(List_lowerBound)  # X轴的数据
    y = np.asarray(List_selectNumber)  # Y轴的数据
    X, Y = np.meshgrid(x, y)  # 生成网格数据
    Z = np.random.rand(len(List_lowerBound), len(List_lowerBound))  # Z轴的数据，这里使用随机数代替
    for i, lowerBound in enumerate(List_lowerBound):
        for j, selectNumber in enumerate(List_selectNumber):
            print(type((lowerBound, selectNumber)))
            print((lowerBound, selectNumber))
            Z[j, i] = CPsRMSET_ParaMining[(lowerBound, selectNumber)]

    # plot figure
    figure_bar(X, Y, Z, path_Figure1)
