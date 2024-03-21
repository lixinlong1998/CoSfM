import Metashape
import csv
import time
import os
import sys

'''
Metashape reconstruction (Only for version 1.8)
Mode: 'FastDSM', 'FastDOM', 'FullyDSM', 'FullyDOM', 'Meshes', 'TiledModel'
'''
#################################################       SETUP      #####################################################
version = 1.8
project_path = ''
chunk_name = ''

AlignImages = True
ErrorReduction = True
BuildDenseCloud = True
BuildDSM = True
WhiteBalance = True
BuildDOM = True
BuildTiledModel = False

ExportDSM = True
ExportDOM = True
ExportTiledModel = True
ExportReport = True
GenerateQuality = True

save_project_each_step = True

# AlignImages Arguments
quality_SPC = 1  # 1(UltraHigh), 2(High), 4(Medium)
Keypoint_Limit = 60000
Keypoint_Limit_Per_Megapixel = 1000
Tiepoint_Limit = 0

# ErrorReduction Arguments
Reprojection_Error = 0.3
Reconstruction_Uncertainty = 15
Projection_Accuracy = 5
MaxRemovePoints = 0.3

# Dense Matching Arguments
quality_DPC = 2  # 1(UltraHigh), 2(High), 4(Medium)
DenseCloud_FileterMode = 'Mild'

# TiledModel Arguments
TiledModel_TiledModelFormat_list = ['Cesium', 'OSGB']
TiledModel_pixel_size = 0.3  # metres
TiledModel_face_count = 20000
TiledModel_ModelFormat = 'None'
TiledModel_ImageFormat = 'JPEG'

# Bundle Adjustment Arguments
bundle_adjustment_args = {
    "fit_f": True,
    "fit_cx": True, "fit_cy": True,
    "fit_b1": False, "fit_b2": False,
    "fit_k1": True, "fit_k2": True,
    "fit_k3": True, "fit_k4": False,
    "fit_p1": True, "fit_p2": True,
    "fit_corrections": False,
    "adaptive_fitting": True,
    "tiepoint_covariance": True
}
#################################################   END OF SETUP   #####################################################
# Checking compatibility
compatible_major_version = str(version)
found_major_version = ".".join(Metashape.app.version.split('.')[:2])
if found_major_version != compatible_major_version:
    raise Exception("Incompatible Metashape version: {} != {}".format(found_major_version, compatible_major_version))
CriteriasThreshold = [Reprojection_Error, Reconstruction_Uncertainty, Projection_Accuracy, MaxRemovePoints]

# some classes of Metashape parameters
TiledModelFormat = {'None': Metashape.TiledModelFormatNone,
                    'TLS': Metashape.TiledModelFormatTLS,
                    'LOD': Metashape.TiledModelFormatLOD,
                    'ZIP': Metashape.TiledModelFormatZIP,
                    'Cesium': Metashape.TiledModelFormatCesium,
                    'SLPK': Metashape.TiledModelFormatSLPK,
                    'OSGB': Metashape.TiledModelFormatOSGB,
                    'OSGT': Metashape.TiledModelFormatOSGT}
ModelFormat = {'None': Metashape.ModelFormatNone,
               'OBJ': Metashape.ModelFormatOBJ,
               '3DS': Metashape.ModelFormat3DS,
               'VRML': Metashape.ModelFormatVRML,
               'PLY': Metashape.ModelFormatPLY,
               'COLLADA': Metashape.ModelFormatCOLLADA,
               'U3D': Metashape.ModelFormatU3D,
               'PDF': Metashape.ModelFormatPDF,
               'DXF': Metashape.ModelFormatDXF,
               'FBX': Metashape.ModelFormatFBX,
               'KMZ': Metashape.ModelFormatKMZ,
               'CTM': Metashape.ModelFormatCTM,
               'STL': Metashape.ModelFormatSTL,
               'DXF_3DF': Metashape.ModelFormatDXF_3DF,
               'TLS': Metashape.ModelFormatTLS,
               'ABC': Metashape.ModelFormatABC,
               'OSGB': Metashape.ModelFormatOSGB,
               'OSGT': Metashape.ModelFormatOSGT,
               'GLTF': Metashape.ModelFormatGLTF,
               'X3D': Metashape.ModelFormatX3D,
               'LandXML': Metashape.ModelFormatLandXML}
ImageFormat = {'None': Metashape.ImageFormatNone,
               'JPEG': Metashape.ImageFormatJPEG,
               'TIFF': Metashape.ImageFormatTIFF,
               'PNG': Metashape.ImageFormatPNG,
               'BMP': Metashape.ImageFormatBMP,
               'EXR': Metashape.ImageFormatEXR,
               'PNM': Metashape.ImageFormatPNM,
               'SGI': Metashape.ImageFormatSGI,
               'CR2': Metashape.ImageFormatCR2,
               'SEQ': Metashape.ImageFormatSEQ,
               'BIL': Metashape.ImageFormatBIL,
               'XYZ': Metashape.ImageFormatXYZ,
               'ARA': Metashape.ImageFormatARA,
               'TGA': Metashape.ImageFormatTGA,
               'DDS': Metashape.ImageFormatDDS,
               'JP2': Metashape.ImageFormatJP2,
               'WebP': Metashape.ImageFormatWebP,
               'JXL': Metashape.ImageFormatJXL}
RasterTransformType = {'None': Metashape.RasterTransformNone,
                       'Value': Metashape.RasterTransformValue,
                       'Palette': Metashape.RasterTransformPalette}
FilterMode = {'No': Metashape.NoFiltering,
              'Mild': Metashape.MildFiltering,
              'Moderate': Metashape.ModerateFiltering,
              'Aggressive': Metashape.AggressiveFiltering}
ImageCompression = {'None': Metashape.ImageCompression.TiffCompressionNone,
                    'LZW': Metashape.ImageCompression.TiffCompressionLZW,
                    'JPEG': Metashape.ImageCompression.TiffCompressionJPEG,
                    'Packbits': Metashape.ImageCompression.TiffCompressionPackbits,
                    'Deflate': Metashape.ImageCompression.TiffCompressionDeflate}

# image compression parameters setting for big tiff exporting
image_compression_paras = Metashape.ImageCompression()
image_compression_paras.tiff_compression = ImageCompression['None']
image_compression_paras.jpeg_quality = 99  # will be only used, if TiffCompressionJPEG is selected#
image_compression_paras.tiff_big = True
image_compression_paras.tiff_overviews = True
image_compression_paras.tiff_tiled = True


def alignImages(chunk, quality=1, KeypointLimit=40000, KeypointLimitPerMegapix=1000, TiepointLimit=4000, realign=True):
    '''
    downscale:quality = 1(UltraHigh), 2(High), 4(Medium)
    '''
    chunk.matchPhotos(downscale=quality, generic_preselection=True, reference_preselection=True,
                      filter_stationary_points=True,
                      keypoint_limit=KeypointLimit, keypoint_limit_per_mpx=KeypointLimitPerMegapix,
                      tiepoint_limit=TiepointLimit, keep_keypoints=True, guided_matching=False, reset_matches=realign)
    chunk.alignCameras(reset_alignment=False)


def bundleAdjustment(chunk, args):
    if version == 1.8:
        chunk.optimizeCameras(fit_f=args['fit_f'],
                              fit_cx=args['fit_cx'], fit_cy=args['fit_cy'],
                              fit_b1=args['fit_b1'], fit_b2=args['fit_b2'],
                              fit_k1=args['fit_k1'], fit_k2=args['fit_k2'],
                              fit_k3=args['fit_k3'], fit_k4=args['fit_k4'],
                              fit_p1=args['fit_p1'], fit_p2=args['fit_p2'],
                              fit_corrections=args['fit_corrections'],
                              adaptive_fitting=args['adaptive_fitting'],
                              tiepoint_covariance=args['tiepoint_covariance'])
    elif version == 1.5:
        chunk.optimizeCameras(fit_f=args['fit_f'],
                              fit_cx=args['fit_cx'], fit_cy=args['fit_cy'],
                              fit_b1=args['fit_b1'], fit_b2=args['fit_b2'],
                              fit_k1=args['fit_k1'], fit_k2=args['fit_k2'],
                              fit_k3=args['fit_k3'], fit_k4=args['fit_k4'],
                              fit_p1=args['fit_p1'], fit_p2=args['fit_p2'],
                              fit_corrections=args['fit_corrections'],
                              adaptive_fitting=args['adaptive_fitting'],
                              tiepoint_covariance=args['tiepoint_covariance'])


def errorReduction(chunk, CriteriasThreshold, bundle_adjustment_args):
    '''
    input:
        Metashape.document.chunk
        CriteriasThreshold =[0.3, 15, 5], refers to RE,RU,PA respectively
    output:
        removedPointIdList=[point_id,...]
    note:
        The removal will not take effect immediately until optimization is performed or the project file is saved
    '''
    print('[Script]    number of points before Error Reduction:', len(chunk.point_cloud.points))
    points = chunk.point_cloud.points
    CriteriasThreshold_RE = CriteriasThreshold[0]
    CriteriasThreshold_RU = CriteriasThreshold[1]
    CriteriasThreshold_PA = CriteriasThreshold[2]
    MaxRemovePoints = CriteriasThreshold[3]

    # Initialize point quality
    # RE, RU, PA, IC
    MF = Metashape.PointCloud.Filter()
    MF.init(chunk, Metashape.PointCloud.Filter.ReprojectionError)
    PointsCriterion_RE = MF.values
    MF.init(chunk, Metashape.PointCloud.Filter.ReconstructionUncertainty)
    PointsCriterion_RU = MF.values
    MF.init(chunk, Metashape.PointCloud.Filter.ProjectionAccuracy)
    PointsCriterion_PA = MF.values
    # MF.init(chunk, Metashape.PointCloud.Filter.ImageCount)
    # PointsCriterion_IC = MF.values

    # PointsQuality
    RemovedPointIdList = []
    for point_id, point in enumerate(points):
        if PointsCriterion_RE[point_id] > CriteriasThreshold_RE:
            RemovedPointIdList.append(point_id)
        if PointsCriterion_RU[point_id] > CriteriasThreshold_RU:
            RemovedPointIdList.append(point_id)
        if PointsCriterion_PA[point_id] > CriteriasThreshold_PA:
            RemovedPointIdList.append(point_id)
    RemovedPointIdSet = set(RemovedPointIdList)
    if len(RemovedPointIdSet) <= len(points) * 0.3:
        # remove outlier points
        for point_id in RemovedPointIdSet:
            point = chunk.point_cloud.points[point_id]
            point.valid = False
    else:
        # in case there are too much points selected to be removed
        while len(RemovedPointIdSet) > len(points) * MaxRemovePoints:
            CriteriasThreshold_RE += 0.1
            CriteriasThreshold_RU += 2
            CriteriasThreshold_PA += 1
            RemovedPointIdList = []
            for point_id, point in enumerate(points):
                if PointsCriterion_RE[point_id] > CriteriasThreshold_RE:
                    RemovedPointIdList.append(point_id)
                if PointsCriterion_RU[point_id] > CriteriasThreshold_RU:
                    RemovedPointIdList.append(point_id)
                if PointsCriterion_PA[point_id] > CriteriasThreshold_PA:
                    RemovedPointIdList.append(point_id)
            RemovedPointIdSet = set(RemovedPointIdList)
        # remove outlier points
        for point_id in RemovedPointIdSet:
            point = chunk.point_cloud.points[point_id]
            point.valid = False
    # Bundle Adjustment,Estimate tiepoint covariance
    bundleAdjustment(chunk, bundle_adjustment_args)
    print('[Script]    number of points after ErrorReduction and optimizeCameras:', len(chunk.point_cloud.points))


def buildDenseCloud(chunk, quality, DenseCloud_FileterMode):
    '''
    downscale:quality = 1(UltraHigh), 2(High), 4(Medium)
    '''
    # build depth maps with moderate filter
    chunk.buildDepthMaps(downscale=quality, filter_mode=FilterMode[DenseCloud_FileterMode])
    # build dense cloud
    chunk.buildDenseCloud(point_colors=True, keep_depth=True, point_confidence=True)
    # export
    # outpath = Metashape.app.document.path[:-4]
    # chunk.exportPoints(path=str(outpath + "_" + str(chunk.label) + "_densecloud.laz"),
    #                    sourceData=Metashape.DataSource.DenseCloudData, save_colors=True, save_confidence=True)


def buildDSM(chunk):
    chunk.buildDem(source_data=Metashape.DataSource.DenseCloudData, interpolation=Metashape.EnabledInterpolation)


def buildDOM(chunk):
    if version == 1.8:
        chunk.buildOrthomosaic(surface_data=Metashape.ElevationData, blending_mode=Metashape.MosaicBlending,
                               fill_holes=True,
                               ghosting_filter=True, cull_faces=False, refine_seamlines=False)
    elif version == 1.5:
        chunk.buildOrthomosaic(surface=Metashape.ElevationData, blending=Metashape.MosaicBlending, fill_holes=True,
                               cull_faces=False, refine_seamlines=False)


def buildTiledModel(chunk, TiledModel_pixel_size, TiledModel_face_count):
    if version == 1.8:
        chunk.buildTiledModel(source_data=Metashape.DataSource.DenseCloudData,
                              pixel_size=TiledModel_pixel_size, tile_size=256, face_count=TiledModel_face_count,
                              ghosting_filter=False, transfer_texture=False, keep_depth=True, merge=False)
    elif version == 1.5:
        chunk.buildTiledModel(tile_size=256, face_count=TiledModel_face_count,
                              ghosting_filter=False, transfer_texture=False, keep_depth=True)


def exportDSM(chunk, path_to_save):
    path_DSM = path_to_save + "_DSM.tif"
    if version == 1.8:
        chunk.exportRaster(path=path_DSM, source_data=Metashape.ElevationData, nodata_value=-32767,
                           save_kml=False, save_world=True, save_scheme=False, save_alpha=False,
                           image_compression=image_compression_paras,
                           image_description=os.path.basename(path_to_save), white_background=False,
                           title='Digital Surface Model')


def exportDOM(chunk, path_to_save):
    path_DOM = path_to_save + "_DOM.tif"
    if version == 1.8:
        chunk.exportRaster(path=path_DOM, source_data=Metashape.OrthomosaicData, nodata_value=-32767,
                           save_kml=False, save_world=True, save_scheme=False, save_alpha=False,
                           image_compression=image_compression_paras,
                           image_description=os.path.basename(path_to_save), white_background=False,
                           title='Orthomosaic')


def exportTiledModel(chunk, path_to_save, TiledModel_TiledModelFormat, TiledModel_ModelFormat, TiledModel_ImageFormat):
    path_TM = path_to_save + f"_TiledModel_{TiledModel_TiledModelFormat}.zip"
    if version == 1.8:
        chunk.exportTiledModel(path=path_TM,
                               format=TiledModelFormat[TiledModel_TiledModelFormat],
                               model_format=ModelFormat[TiledModel_ModelFormat],
                               texture_format=ImageFormat[TiledModel_ImageFormat],
                               raster_transform=Metashape.RasterTransformNone,
                               clip_to_boundary=True, model_compression=True, screen_space_error=16)


def exportReport(chunk, path_to_save):
    path_report = path_to_save + '_report.pdf'
    name_document = os.path.basename(path_to_save)
    chunk.exportReport(path=path_report, title=name_document)


def exportQuality(chunk, path_to_save):
    path_quality = path_to_save + '_quality.txt'
    f = open(path_quality, "w")
    fwriter = csv.writer(f, delimiter='\t', lineterminator='\n')
    fwriter.writerow(str(chunk.orthomosaic.resolution))
    fwriter.writerow(str(chunk.elevation.resolution))
    f.close()


starttime0 = time.perf_counter()
if __name__ == '__main__':
    # [1]  open document
    if project_path:
        # run script from cmd
        doc = Metashape.app.document
        doc.open(project_path)
    else:
        # run script from GUI
        doc = Metashape.app.document
    project_path = doc.path
    workspace = os.path.dirname(project_path)
    project_name = os.path.basename(project_path)
    project_name_no_end = os.path.splitext(project_name)[0]

    # [2]  access chunks
    chunksList = []
    if chunk_name:
        # choose the chunk with given name
        for chunk_i in Metashape.app.document.chunks:
            if chunk_i.label == chunk_name:
                chunksList.append(chunk_i)
    else:
        # choose the chunk enabled
        for chunk_i in Metashape.app.document.chunks:
            if chunk_i.enabled:
                chunksList.append(chunk_i)
            else:
                continue

    # [3]  access chunk label
    for chunk_id, chunk in enumerate(chunksList):
        print(chunk.label)
        path_to_save = workspace + f"/{project_name_no_end}_{chunk_id}"

        # [4] build alignCameras()
        if AlignImages and not chunk.point_cloud:
            alignImages(chunk, quality=quality_SPC, KeypointLimit=Keypoint_Limit,
                        KeypointLimitPerMegapix=Keypoint_Limit_Per_Megapixel, TiepointLimit=Tiepoint_Limit)
            if save_project_each_step:
                doc.save()

        # [5] gradual remove outliers
        if ErrorReduction and not chunk.dense_cloud:
            errorReduction(chunk, CriteriasThreshold, bundle_adjustment_args)
            if save_project_each_step:
                doc.save()

        # [6]  build dense cloud for each individual chunks
        if BuildDenseCloud and not chunk.dense_cloud:
            buildDenseCloud(chunk, quality_DPC, DenseCloud_FileterMode)
            if save_project_each_step:
                doc.save()

        # [7]  build DSM
        if BuildDSM and not chunk.elevation:
            buildDSM(chunk)
            if save_project_each_step:
                doc.save()

        # [8]  calibrate colors and white balance
        if WhiteBalance and chunk.elevation and not chunk.orthomosaic:
            chunk.calibrateColors(source_data=Metashape.DataSource.ElevationData, white_balance=True)
            if save_project_each_step:
                doc.save()

        # [9] build DOM
        if BuildDOM and chunk.elevation and not chunk.orthomosaic:
            buildDOM(chunk)
            if save_project_each_step:
                doc.save()

        # [10] build TiledModel
        if BuildTiledModel and chunk.elevation and not chunk.tiled_model:
            buildTiledModel(chunk, TiledModel_pixel_size, TiledModel_face_count)
            if save_project_each_step:
                doc.save()

        # [11] export DSM, DOM and TiledModel
        if chunk.elevation and ExportDSM:
            exportDSM(chunk, path_to_save)
        if chunk.orthomosaic and ExportDOM:
            exportDOM(chunk, path_to_save)
        if chunk.tiled_model and ExportTiledModel:
            for TiledModel_TiledModelFormat in TiledModel_TiledModelFormat_list:
                exportTiledModel(chunk, path_to_save, TiledModel_TiledModelFormat, TiledModel_ModelFormat,
                                 TiledModel_ImageFormat)
        # [12]  export report
        if ExportReport:
            exportReport(chunk, path_to_save)
        # [13]
        if GenerateQuality:
            exportQuality(chunk, path_to_save)
        # save project
        doc.save()
    ProcessTime = time.perf_counter() - starttime0
    print('[Script]    All processed done in {0}h {1}m {2:0.3f} sec.'.format(ProcessTime // 3600,
                                                                         ProcessTime % 3600 // 60,
                                                                         ProcessTime % 60))
