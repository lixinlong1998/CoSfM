# coding:utf-8
import os
import sys
import time

'''
Before you start using this plugin, please do:
    1.set up the following paths.
        PATH_CODE = 'the root path of this code'
        PATH_CUDA = 'the location of bin folder of your cuda'
        metashape_version = this code only support 1.8
    2.Make sure the library functions imported below are installed.
    3.If you want to use the SIFT with cuda, please set up the 'PATH_CODE' and 'PATH_OPENCV'
'''
########################################################################################################################
################################################       SETUP       #####################################################
metashape_version = 1.8
PATH_CODE = r'D:\Research\20221223_CoSfM\Release\CFTM_v1.0'
PATH_COLMAP_BAT = r'D:\Software\COLMAP\COLMAP-3.6-windows-cuda\COLMAP.bat'
# PATH_CUDA = r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8\bin'
# PATH_OPENCV = r'G:\UAV_SOFTWARE\OpenCV\opencv_4_5_0_cuda_11_1_py38\install\x64\vc16\bin'
###################################################   END SETUP   ######################################################
########################################################################################################################

PATH_SRCCOSFM = os.path.join(PATH_CODE, 'src_CFTM')
PATH_SRCMETASHAPE = os.path.join(PATH_CODE, 'src_Metashape')
PATH_TOOLBOX = os.path.join(PATH_CODE, 'src_Process')
sys.path.append(PATH_CODE)
sys.path.append(PATH_SRCCOSFM)
sys.path.append(PATH_SRCMETASHAPE)
sys.path.append(PATH_TOOLBOX)
# os.add_dll_directory(PATH_CUDA)
# os.add_dll_directory(PATH_OPENCV)

# import math
# import sqlite3
# import random
# import multiprocessing
# import scipy
# import numpy as np
# import cv2 as cv
# import matplotlib.pyplot as plt
# import networkx as nx
# from osgeo import ogr, osr
# from osgeo import gdal
# from scipy.stats import gaussian_kde
# from scipy.optimize import least_squares
