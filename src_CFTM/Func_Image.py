import os
import time
import sys
import csv
import random
import numpy as np
import matplotlib.pyplot as plt
from osgeo import ogr, osr
from osgeo import gdal
import numpy as np
import src_CFTM.Func_Geometry as FcGeo


def readImg(path):
    # 读取栅格数据
    tif = gdal.Open(path)
    # 获取栅格文件的地理转换和波段
    geo_transform = tif.GetGeoTransform()
    projection = tif.GetProjection()
    width = tif.RasterXSize  # 栅格矩阵的列数
    height = tif.RasterYSize  # 栅格矩阵的行数
    bands = tif.RasterCount  # 波段数
    # tif_b1 = tif.GetRasterBand(1)
    # 将tif转换为数组
    tif_array = tif.ReadAsArray()
    return tif_array


def writeImg(path, im_proj, im_geotrans, data_array):
    if 'int8' in data_array.dtype.name:
        datatype = gdal.GDT_Int16
    elif 'int16' in data_array.dtype.name:
        datatype = gdal.GDT_Int16
    else:
        datatype = gdal.GDT_Float32

    if len(data_array.shape) == 3:
        im_bands, im_height, im_width = data_array.shape
    else:
        im_bands, (im_height, im_width) = 1, data_array.shape

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_height, im_bands, datatype)

    dataset.SetGeoTransform(im_geotrans)
    dataset.SetProjection(im_proj)

    if im_bands == 1:
        dataset.GetRasterBand(1).WriteArray(data_array)
    else:
        for i in range(im_bands):
            dataset.GetRasterBand(i + 1).WriteArray(data_array[i])

    del dataset


def maskCTPs_Unstable(MarkersGrid, unstable_areas_mask_path):
    # UnstableRegion = ogr.Open(unstable_areas_mask_path)
    driver = ogr.GetDriverByName("ESRI Shapefile")
    # <osgeo.ogr.DataSource; proxy of <Swig Object of type 'OGRDataSourceShadow *' at 0x000001F7359F7B70> >
    UnstableRegion = driver.Open(unstable_areas_mask_path, 1)
    MarkersNominated = {}
    for grid_id, Markers in MarkersGrid.items():
        MarkersNominated_grid = []
        for i, Marker in enumerate(Markers):
            Point_Coord_PJCS_XY = [Marker[0][1][0], Marker[0][1][1]]
            if FcGeo.isInPolygonShp(Point_Coord_PJCS_XY, UnstableRegion):
                continue
            MarkersNominated_grid.append(Marker)
        MarkersNominated[grid_id] = MarkersNominated_grid
    return MarkersNominated


def plot_hist_together(path, data, bin_width=0.5):
    flattened1 = data[0].flatten()
    flattened2 = data[1].flatten()
    flattened3 = data[2].flatten()
    data1 = flattened1[~np.isnan(flattened1)]
    data2 = flattened2[~np.isnan(flattened2)]
    data3 = flattened3[~np.isnan(flattened3)]
    print(data1)
    # 计算每个指标的最小值、最大值和范围,根据给定的柱子宽度计算柱子的边界值
    x_bins = np.arange(np.min(data1), np.max(data1) + bin_width, bin_width)
    y_bins = np.arange(np.min(data2), np.max(data2) + bin_width, bin_width)
    z_bins = np.arange(np.min(data3), np.max(data3) + bin_width, bin_width)
    # 创建一个大小为6x4英寸的画布和一个Axes对象
    fig, ax = plt.subplots(figsize=(10, 6), dpi=128)
    # 绘制直方图
    ax.hist(data1, bins=x_bins, alpha=0.3, color='blue', label='Error(pixel)')
    ax.hist(data2, bins=y_bins, alpha=0.3, color='green', label='RepErrV3e1')
    ax.hist(data3, bins=z_bins, alpha=0.3, color='red', label='RepErrV3e2')
    # 添加图例和坐标轴标签
    ax.legend(loc='best')
    ax.set_xlabel('index')
    ax.set_ylabel('frequency')
    # 调整子图之间的距离和周围留白
    plt.subplots_adjust(wspace=0.3, left=0.1, right=0.95, bottom=0.15)
    # 显示图形
    plt.savefig(path)
    plt.close()


def clip_raster_with_shapefile(input_shp, input_tif):
    # 打开矢量文件和栅格文件
    shp_ds = ogr.Open(input_shp)
    tif_ds = gdal.Open(input_tif)
    output_tif = input_tif.replace('.tif', '_Stable.tif')

    # 获取矢量文件的几何信息和投影
    layer = shp_ds.GetLayer()
    feature = layer.GetNextFeature()
    geometry = feature.GetGeometryRef()
    spatial_ref = layer.GetSpatialRef()

    # 获取栅格文件的地理转换和波段
    geo_transform = tif_ds.GetGeoTransform()
    projection = tif_ds.GetProjection()
    band = tif_ds.GetRasterBand(1)

    # 创建输出栅格文件
    driver = gdal.GetDriverByName('GTiff')
    output_ds = driver.Create(output_tif, band.XSize, band.YSize, 1, band.DataType)
    output_ds.SetGeoTransform(geo_transform)
    output_ds.SetProjection(projection)

    # 将矢量文件的几何信息转换为栅格文件的像素坐标系
    transform = osr.CoordinateTransformation(spatial_ref, osr.SpatialReference(wkt=output_ds.GetProjection()))
    geometry.Transform(transform)

    # 通过裁剪栅格文件中的像素来创建新的栅格图像
    output_band = output_ds.GetRasterBand(1)
    output_band.SetNoDataValue(0)  # 设置无效像素值
    output_band.FlushCache()

    gdal.RasterizeLayer(output_ds, [1], layer, burn_values=[1])

    # 将原始影像的地理转换和投影信息应用于结果影像
    output_ds.SetGeoTransform(geo_transform)
    output_ds.SetProjection(projection)

    # 读取原始影像的像素值并应用于结果影像
    original_data = band.ReadAsArray()
    output_data = output_band.ReadAsArray()
    output_data[output_data == 1] = original_data[output_data == 1]
    output_band.WriteArray(output_data)

    # 将值为0的像素设置为NaN
    output_data[output_data == 0] = np.nan
    # 计算统计指标
    mean = np.nanmean(output_data)
    std = np.nanstd(output_data)
    abs_mean = np.nanmean(np.abs(output_data))
    print(mean, std, abs_mean)
    # 关闭数据集
    shp_ds = None
    tif_ds = None
    return output_data