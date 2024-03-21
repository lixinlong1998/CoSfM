# coding:utf-8
import os
import sys
import time

'''
Before you start using this plugin, please do:
    1.set up the following paths.
        PATH_CODE = 'the root path of this code'
        PATH_CUDA = 'the location of bin folder of your cuda'
        PATH_OPENCV = 'the location of bin folder of your opencv'
        metashape_version = this version only support 1.8.5
    2.Make sure the library functions imported below are installed.
    
'''

metashape_version = 1.8
PATH_CODE = r'G:\AResearchG\20221223_CoSfM\Release\CFTM_v1.0'
PATH_CUDA = r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8\bin'
PATH_OPENCV = r'G:\UAV_SOFTWARE\OpenCV\opencv_4_5_0_cuda_11_1_py38\install\x64\vc16\bin'
PATH_COLMAP_BAT = r'G:\UAV_SOFTWARE\Colmap\COLMAP-3.6-windows-cuda_PreRelease\COLMAP-3.6-windows-cuda_2\colmap.bat'

PATH_SRCCOSFM = os.path.join(PATH_CODE, 'src_CFTM')
PATH_SRCMETASHAPE = os.path.join(PATH_CODE, 'src_Metashape')
PATH_TOOLBOX = os.path.join(PATH_CODE, 'src_Process')

sys.path.append(PATH_CODE)
sys.path.append(PATH_SRCCOSFM)
sys.path.append(PATH_SRCMETASHAPE)
sys.path.append(PATH_TOOLBOX)
os.add_dll_directory(PATH_CUDA)
os.add_dll_directory(PATH_OPENCV)

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
