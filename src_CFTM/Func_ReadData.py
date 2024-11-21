# coding:gbk
import numpy as np
import csv
import os
import sys


def readMarkersAllQua(path):
    '''
    input:
        Qua.txt
    output:
        Qua = [i,
               0
               Point3D_e1_PJCS[0], Point3D_e1_PJCS[1], Point3D_e1_PJCS[2],
               1                   2                   3
               Point3D_e2_PJCS[0], Point3D_e2_PJCS[1], Point3D_e2_PJCS[2],
               4                   5                   6
               residual_X, residual_Y, residual_Z, residual_P, residual_T,
               7           8           9           10          11
               error_pixel, RepErr_v3_e1, RepErr_v3_e2, RepErr_v2_e1, RepErr_v2_e2,
               12           13            14            15            16
               TriAng_Max, TriAng_Avg, TriAng_Min,
               17          18          19
               ViewNum, ViewNum1, ViewNum2, EpoMatch_Num]
               20       21        22        23
    '''
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        head = next(freader)  # skip header
        MarkersQua = []
        for row in freader:
            MarkersQua.append([eval(item) for item in row])
    return MarkersQua


def readMarkersAllQuaReport(path):
    '''
    input:
        QuaSta.txt
    output:
        QuaSta={}
    Res_X_RMSE(m)   Res_Y_RMSE(m)   Res_Z_RMSE(m)   Res_P_RMSE(m)   Res_T_RMSE(m)
    0               1               2               3               4
    Res_X_MAE(m)    Res_Y_MAE(m)    Res_Z_MAE(m)    Res_P_MAE(m)    Res_T_MAE(m)
    5               6               7               8               9
    Res_X_AVG(m)    Res_Y_AVG(m)    Res_Z_AVG(m)    Res_P_AVG(m)    Res_T_AVG(m)
    10              11              12              13              14
    Res_X_STD(m)    Res_Y_STD(m)    Res_Z_STD(m)    Res_P_STD(m)    Res_T_STD(m)
    15              16              17              18              19
    Distrib_residual_X_AVG      Distrib_residual_X_MIN      Distrib_residual_X_25P      Distrib_residual_X_50P      Distrib_residual_X_75P      Distrib_residual_X_MAX
    20                          21                          22                          23                          24                          25
    Distrib_residual_Z_AVG      Distrib_residual_Z_MIN      Distrib_residual_Z_25P      Distrib_residual_Z_50P      Distrib_residual_Z_75P      Distrib_residual_Z_MAX
    26
    Distrib_residual_Y_AVG      Distrib_residual_Y_MIN      Distrib_residual_Y_25P      Distrib_residual_Y_50P      Distrib_residual_Y_75P      Distrib_residual_Y_MAX
    32
    Distrib_residual_P_AVG      Distrib_residual_P_MIN      Distrib_residual_P_25P      Distrib_residual_P_50P      Distrib_residual_P_75P      Distrib_residual_P_MAX
    38
    Distrib_residual_T_AVG      Distrib_residual_T_MIN      Distrib_residual_T_25P      Distrib_residual_T_50P      Distrib_residual_T_75P      Distrib_residual_T_MAX
    44
    Distrib_error_pixel_AVG     Distrib_error_pixel_MIN     Distrib_error_pixel_25P     Distrib_error_pixel_50P     Distrib_error_pixel_75P     Distrib_error_pixel_MAX
    50
    Distrib_RepErr_v3_e1_AVG    Distrib_RepErr_v3_e1_MIN    Distrib_RepErr_v3_e1_25P    Distrib_RepErr_v3_e1_50P    Distrib_RepErr_v3_e1_75P    Distrib_RepErr_v3_e1_MAX
    56
    Distrib_RepErr_v3_e2_AVG    Distrib_RepErr_v3_e2_MIN    Distrib_RepErr_v3_e2_25P    Distrib_RepErr_v3_e2_50P    Distrib_RepErr_v3_e2_75P    Distrib_RepErr_v3_e2_MAX
    62
    Distrib_RepErr_v2_e1_AVG    Distrib_RepErr_v2_e1_MIN    Distrib_RepErr_v2_e1_25P    Distrib_RepErr_v2_e1_50P    Distrib_RepErr_v2_e1_75P    Distrib_RepErr_v2_e1_MAX
    68
    Distrib_RepErr_v2_e2_AVG    Distrib_RepErr_v2_e2_MIN    Distrib_RepErr_v2_e2_25P    Distrib_RepErr_v2_e2_50P    Distrib_RepErr_v2_e2_75P    Distrib_RepErr_v2_e2_MAX
    74
    Distrib_TriAng_Max_AVG      Distrib_TriAng_Max_MIN      Distrib_TriAng_Max_25P      Distrib_TriAng_Max_50P      Distrib_TriAng_Max_75P      Distrib_TriAng_Max_MAX
    80
    Distrib_TriAng_Avg_AVG      Distrib_TriAng_Avg_MIN      Distrib_TriAng_Avg_25P      Distrib_TriAng_Avg_50P      Distrib_TriAng_Avg_75P      Distrib_TriAng_Avg_MAX
    86
    Distrib_TriAng_Min_AVG      Distrib_TriAng_Min_MIN      Distrib_TriAng_Min_25P      Distrib_TriAng_Min_50P      Distrib_TriAng_Min_75P      Distrib_TriAng_Min_MAX
    92
    Distrib_ViewNum_AVG         Distrib_ViewNum_MIN         Distrib_ViewNum_25P         Distrib_ViewNum_50P         Distrib_ViewNum_75P         Distrib_ViewNum_MAX
    98
    Distrib_ViewNum1_AVG        Distrib_ViewNum1_MIN        Distrib_ViewNum1_25P        Distrib_ViewNum1_50P        Distrib_ViewNum1_75P        Distrib_ViewNum1_MAX
    104
    Distrib_ViewNum2_AVG        Distrib_ViewNum2_MIN        Distrib_ViewNum2_25P        Distrib_ViewNum2_50P        Distrib_ViewNum2_75P        Distrib_ViewNum2_MAX
    110
    Distrib_EpoMatch_Num_AVG    Distrib_EpoMatch_Num_MIN    Distrib_EpoMatch_Num_25P    Distrib_EpoMatch_Num_50P    Distrib_EpoMatch_Num_75P    Distrib_EpoMatch_Num_MAX
    116                         117                         118                         119                         120                         121
    '''
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        fileData = []
        MarkersQuaSta = {}
        for row in freader:
            fileData.append(row)
        for i in range(0, 8, 2):
            for j in range(len(fileData[i])):
                MarkersQuaSta[fileData[i][j]] = eval(fileData[i + 1][j])
        for i in range(8, len(fileData), 2):
            MarkersQuaSta[fileData[i][0][0:-4] + '_AVG'] = eval(fileData[i + 1][0])
            MarkersQuaSta[fileData[i][0][0:-4] + '_MIN'] = eval(fileData[i + 1][1])
            MarkersQuaSta[fileData[i][0][0:-4] + '_25P'] = eval(fileData[i + 1][2])
            MarkersQuaSta[fileData[i][0][0:-4] + '_50P'] = eval(fileData[i + 1][3])
            MarkersQuaSta[fileData[i][0][0:-4] + '_75P'] = eval(fileData[i + 1][4])
            MarkersQuaSta[fileData[i][0][0:-4] + '_MAX'] = eval(fileData[i + 1][5])
    return MarkersQuaSta


def mergeCPsCTPsQuaSta(path_RandSelect, subFolder_RandomSelect):
    '''
    input:
        QuaSta.txt
    output:
        Merged QuaSta of CPs and CTPs, where each index contain two values from CPsQuaSta and CTPsQuaSta respectively.
    Res_X_RMSE(m)   Res_Y_RMSE(m)   Res_Z_RMSE(m)   Res_P_RMSE(m)   Res_T_RMSE(m)
        0,1             2,3             4,5             6,7             8,9
    Res_X_MAE(m)    Res_Y_MAE(m)    Res_Z_MAE(m)    Res_P_MAE(m)    Res_T_MAE(m)
       10,11           12,13           14,15           16,17           18,19
    Res_X_AVG(m)    Res_Y_AVG(m)    Res_Z_AVG(m)    Res_P_AVG(m)    Res_T_AVG(m)
       20,21           22,23           24,25           26,27           28,29
    Res_X_STD(m)    Res_Y_STD(m)    Res_Z_STD(m)    Res_P_STD(m)    Res_T_STD(m)
       30,31           32,33           34,35           36,37           38,39
    Distrib_residual_X_AVG      Distrib_residual_X_MIN      Distrib_residual_X_25P      Distrib_residual_X_50P      Distrib_residual_X_75P      Distrib_residual_X_MAX
    40,41                       42,43                       44,45                       46,47                       48,49                       50,51
    Distrib_residual_Z_AVG      Distrib_residual_Z_MIN      Distrib_residual_Z_25P      Distrib_residual_Z_50P      Distrib_residual_Z_75P      Distrib_residual_Z_MAX
    52,53                       42+12x1.43+12x1
    Distrib_residual_Y_AVG      Distrib_residual_Y_MIN      Distrib_residual_Y_25P      Distrib_residual_Y_50P      Distrib_residual_Y_75P      Distrib_residual_Y_MAX
    64,65                       42+12x2.43+12x2
    Distrib_residual_P_AVG      Distrib_residual_P_MIN      Distrib_residual_P_25P      Distrib_residual_P_50P      Distrib_residual_P_75P      Distrib_residual_P_MAX
    76,77
    Distrib_residual_T_AVG      Distrib_residual_T_MIN      Distrib_residual_T_25P      Distrib_residual_T_50P      Distrib_residual_T_75P      Distrib_residual_T_MAX
    88,89
    Distrib_Error_pixel_AVG     Distrib_error_pixel_MIN     Distrib_error_pixel_25P     Distrib_error_pixel_50P     Distrib_error_pixel_75P     Distrib_error_pixel_MAX
    100,101
    Distrib_RepErr_v3_e1_AVG    Distrib_RepErr_v3_e1_MIN    Distrib_RepErr_v3_e1_25P    Distrib_RepErr_v3_e1_50P    Distrib_RepErr_v3_e1_75P    Distrib_RepErr_v3_e1_MAX
    112,113
    Distrib_RepErr_v3_e2_AVG    Distrib_RepErr_v3_e2_MIN    Distrib_RepErr_v3_e2_25P    Distrib_RepErr_v3_e2_50P    Distrib_RepErr_v3_e2_75P    Distrib_RepErr_v3_e2_MAX
    124,125
    Distrib_RepErr_v2_e1_AVG    Distrib_RepErr_v2_e1_MIN    Distrib_RepErr_v2_e1_25P    Distrib_RepErr_v2_e1_50P    Distrib_RepErr_v2_e1_75P    Distrib_RepErr_v2_e1_MAX
    136,137
    Distrib_RepErr_v2_e2_AVG    Distrib_RepErr_v2_e2_MIN    Distrib_RepErr_v2_e2_25P    Distrib_RepErr_v2_e2_50P    Distrib_RepErr_v2_e2_75P    Distrib_RepErr_v2_e2_MAX
    148,149
    Distrib_TriAng_Max_AVG      Distrib_TriAng_Max_MIN      Distrib_TriAng_Max_25P      Distrib_TriAng_Max_50P      Distrib_TriAng_Max_75P      Distrib_TriAng_Max_MAX
    160,161
    Distrib_TriAng_Avg_AVG      Distrib_TriAng_Avg_MIN      Distrib_TriAng_Avg_25P      Distrib_TriAng_Avg_50P      Distrib_TriAng_Avg_75P      Distrib_TriAng_Avg_MAX
    172,173
    Distrib_TriAng_Min_AVG      Distrib_TriAng_Min_MIN      Distrib_TriAng_Min_25P      Distrib_TriAng_Min_50P      Distrib_TriAng_Min_75P      Distrib_TriAng_Min_MAX
    184,185
    Distrib_ViewNum_AVG         Distrib_ViewNum_MIN         Distrib_ViewNum_25P         Distrib_ViewNum_50P         Distrib_ViewNum_75P         Distrib_ViewNum_MAX
    196,197
    Distrib_ViewNum1_AVG        Distrib_ViewNum1_MIN        Distrib_ViewNum1_25P        Distrib_ViewNum1_50P        Distrib_ViewNum1_75P        Distrib_ViewNum1_MAX
    208,209
    Distrib_ViewNum2_AVG        Distrib_ViewNum2_MIN        Distrib_ViewNum2_25P        Distrib_ViewNum2_50P        Distrib_ViewNum2_75P        Distrib_ViewNum2_MAX
    220,221
    Distrib_EpoMatch_Num_AVG    Distrib_EpoMatch_Num_MIN    Distrib_EpoMatch_Num_25P    Distrib_EpoMatch_Num_50P    Distrib_EpoMatch_Num_75P    Distrib_EpoMatch_Num_MAX
    232,233                     234,235                     236,237                     238,239                     240,241                     242,243
    '''
    # read file and convert data format
    CPsCTPsQuaSta_Database = []
    CPsCTPsQuaSta_Database_ids = []
    for a, subFolderName in enumerate(subFolder_RandomSelect):
        for root, dirs, files in os.walk(path_RandSelect + '/' + subFolderName):
            print(root)
            # 每一个子文件夹内部的文件序号是从0开始的，所以需要单独对每个子文件夹处理
            pathList_CPsQuaReport = {}
            pathList_CTPsQuaReport = {}
            for fileName in files:
                fileTypeName = fileName.split('_')[0]
                b = int(fileName[0:-4].split('_')[1])  # 文件顺序是按照系统的命名排序来的，即1过了是10、11、而不是2
                if fileTypeName == 'CPsQuaSta':
                    path_CPsQuaSta = root + '/' + fileName
                    pathList_CPsQuaReport[b] = path_CPsQuaSta
                if fileTypeName == 'CTPsQuaSta':
                    path_CTPsQuaSta = root + '/' + fileName
                    pathList_CTPsQuaReport[b] = path_CTPsQuaSta
            # 通过序号找到一次迭代的结果CPs和CTPs，并将它们合并
            for b, path_CPsQuaSta in pathList_CPsQuaReport.items():
                path_CTPsQuaSta = pathList_CTPsQuaReport[b]
                # 根据文件地址读取数据
                CPsQuaReport = readMarkersAllQuaReport(path_CPsQuaSta)
                CTPsQuaReport = readMarkersAllQuaReport(path_CTPsQuaSta)
                # 合并CTPsQuaReport和CPsQuaReport
                CPsCTPsQuaSta = []
                CPsCTPsQuaSta_id = b
                for indexName, CPsQuaValue in CPsQuaReport.items():
                    CTPsQuaValue = CTPsQuaReport[indexName]
                    CPsCTPsQuaSta.append(CPsQuaValue)
                    CPsCTPsQuaSta.append(CTPsQuaValue)
                CPsCTPsQuaSta_Database.append(CPsCTPsQuaSta)
                CPsCTPsQuaSta_Database_ids.append(CPsCTPsQuaSta_id)
    return CPsCTPsQuaSta_Database, CPsCTPsQuaSta_Database_ids


def statistic_analysis(matrix: np.ndarray, i):
    stat_list = []
    stat_list.append(np.mean(matrix, axis=0).tolist())
    stat_list.append(np.std(matrix, axis=0).tolist())
    stat_list.append(np.max(matrix, axis=0).tolist())
    stat_list.append(np.percentile(matrix, 75, axis=0).tolist())
    stat_list.append(np.percentile(matrix, 50, axis=0).tolist())
    stat_list.append(np.percentile(matrix, 25, axis=0).tolist())
    stat_list.append(np.min(matrix, axis=0).tolist())
    stat_matrix = np.mat(stat_list)
    return np.around(stat_matrix, i)


def readQuas(path_OptimalSelct):
    # read file and convert data format
    for root, dirs, files in os.walk(path_OptimalSelct):
        # 每一个子文件夹内部的文件序号是从0开始的，所以需要单独对每个子文件夹处理
        CTPsQua_Dict = {}
        CTPsQuaSta_Dict = {}
        CPsQuaSta_Dict = {}
        for fileName in files:
            fileType = fileName.split('_')[0]
            fileId = int(fileName[0:-4].split('_')[1])
            if fileType == 'CPsQuaSta':
                path = root + '/' + fileName
                CPsQuaSta = readMarkersAllQuaReport(path)
                CPsQuaSta_Dict[fileId] = CPsQuaSta
            if fileType == 'CTPsQua':
                path = root + '/' + fileName
                CTPsQua = readMarkersAllQua(path)
                CTPsQua_Dict[fileId] = CTPsQua
            if fileType == 'CTPsQuaSta':
                path = root + '/' + fileName
                CTPsQuaSta = readMarkersAllQuaReport(path)
                CTPsQuaSta_Dict[fileId] = CTPsQuaSta
    # sort the dictionary to list
    CTPsQua_List = []
    CTPsQuaSta_List = []
    CPsQuaSta_List = []
    for i in range(len(CPsQuaSta_Dict)):
        CTPsQua_List.append(CTPsQua_Dict[i])
        CTPsQuaSta_List.append(CTPsQuaSta_Dict[i])
        CPsQuaSta_List.append(CPsQuaSta_Dict[i])
    return CTPsQua_List, CTPsQuaSta_List, CPsQuaSta_List


def getQuaIndList(CTPsQua_List, iteratNum, Ind=0):
    CTPsQuaIndex_List = []
    for i, CTPsQua in enumerate(CTPsQua_List):
        if i < iteratNum:
            CTPsQuaIndex_List.append([CTPQua[Ind] for CTPQua in CTPsQua])
    return CTPsQuaIndex_List
