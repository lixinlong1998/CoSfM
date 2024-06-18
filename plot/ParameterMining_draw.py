import csv
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.cm import get_cmap
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(r'D:\Research\20221223_CoSfM\Release\CFTM_v1.0')
from src_CFTM import ConnectData_Metashape
from src_Metashape import FuncMs_Basic, FuncMs_Marker
from src_CFTM import Func_ReadData, Func_CommonTiePoints

'''
Please run this script using the following command in cmd：

metashape.exe -r D:\Research\20221223_CoSfM\Release\CFTM_v1.0\plot\ParameterMining_draw.py

'''

#################################################       SETUP      #####################################################
path_OptmSelect = r'E:\Projects\20230418_CFTM\20240509_Revision\Exp5_ParameterMining'
path_Figure1 = r'E:\Projects\20230418_CFTM\20240509_Revision\ParamMining.png'
List_lowerBound = [10, 20, 30, 40, 50, 60]
List_selectNumber = [1, 5, 10, 15, 20, 25]


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
    pathList_CPsQuaReport = {}
    CPsCTPsQuaSta_Database = []
    CPsCTPsQuaSta_Database_ids = []
    for root, dirs, files in os.walk(path_OptmSelect):
        print(files)
        for fileName in files:
            fileTypeName = fileName.split('_')[1]
            grid_id = eval(fileName.split('_')[0])  # 文件顺序是按照系统的命名排序来的，即1过了是10、11、而不是2
            print(grid_id)
            if fileTypeName == 'CPsQuaSta.txt':
                path_CPsQuaSta = root + '/' + fileName
                CPsQuaReport = Func_ReadData.readMarkersAllQuaReport(path_CPsQuaSta)
                CPsCTPsQuaSta = []
                for indexName, CPsQuaValue in CPsQuaReport.items():
                    CPsCTPsQuaSta.append(CPsQuaValue)
                CPsCTPsQuaSta_Database.append(CPsCTPsQuaSta)
                CPsCTPsQuaSta_Database_ids.append(grid_id)

    rank_col_id = 4
    # print(CPsCTPsQuaSta_Database)
    db = np.asarray(CPsCTPsQuaSta_Database)
    db_ids = CPsCTPsQuaSta_Database_ids
    CPsRMSET = list(db[:, rank_col_id])
    print('iterations num:', len(CPsRMSET))
    # CPsRMSE-T最大的CTPs
    print('iterations: Max CPsRMSET:{0} in {1}'.format(max(CPsRMSET), db_ids[CPsRMSET.index(max(CPsRMSET))]))
    # CPsRMSE-T最小的CTPs
    print('iterations: Min CPsRMSET:{0} in {1}'.format(min(CPsRMSET), db_ids[CPsRMSET.index(min(CPsRMSET))]))

    # CPsRMSET_PM_ListK = [CPsRMSET for test_id, CPsRMSET in CPsRMSET_ParaMining.items()]
    # CPsRMSET_PM_ListV = [test_id for test_id, CPsRMSET in CPsRMSET_ParaMining.items()]
    # print(max(CPsRMSET_PM_ListK))
    # print(CPsRMSET_PM_ListV[CPsRMSET_PM_ListK.index(max(CPsRMSET_PM_ListK))])
    # print(min(CPsRMSET_PM_ListK))
    # print(CPsRMSET_PM_ListV[CPsRMSET_PM_ListK.index(min(CPsRMSET_PM_ListK))])
    #
    # # 定义数据
    # x = np.asarray(List_lowerBound)  # X轴的数据
    # y = np.asarray(List_selectNumber)  # Y轴的数据
    # X, Y = np.meshgrid(x, y)  # 生成网格数据
    # Z = np.random.rand(len(List_lowerBound), len(List_lowerBound))  # Z轴的数据，这里使用随机数代替
    # for i, lowerBound in enumerate(List_lowerBound):
    #     for j, selectNumber in enumerate(List_selectNumber):
    #         Z[j, i] = CPsRMSET_ParaMining[str((lowerBound, selectNumber))]
    #
    # # plot figure
    # figure_bar(X, Y, Z, path_Figure1)
