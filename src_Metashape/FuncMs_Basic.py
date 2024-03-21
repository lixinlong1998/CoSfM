from config import *
import Metashape
import csv

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


def openDocument(args):
    project_path = args["project_path"]
    if project_path:
        # open exist document
        document = Metashape.app.document
        document.open(project_path)
    else:
        # creat a new document
        document = Metashape.app.document
        document.save(project_path)

        # add photos
        images_path = args["images_path"]
        images_list = []
        for root, dirs, files in os.walk(images_path):
            for image_folder_path in dirs:
                images_list += [entry.path for entry in os.scandir(image_folder_path) if
                                (entry.is_file() and os.path.splitext(entry.name)[1].lower() in
                                 [".jpg", ".jpeg", ".tif", ".tiff"])]
        chunk = document.addChunk()
        chunk.addPhotos(images_list)

        # set coordinate system
        PJCS_EPSG = args["reference_args"]["PJCS_EPSG"]
        crs_wkt = chunk.crs.wkt  # Coordinate wkt format
        if crs_wkt[:6] != 'PROJCS':
            new_crs = Metashape.CoordinateSystem(f'EPSG::{PJCS_EPSG}')
            chunk.crs = new_crs

        # set accuracy of reference
        camera_accuracy_m = args["reference_args"]["camera_accuracy_m"]
        camera_accuracy_deg = args["reference_args"]["camera_accuracy_deg"]
        marker_accuracy_m = args["reference_args"]["marker_accuracy_m"]
        marker_accuracy_pix = args["reference_args"]["marker_accuracy_pix"]
        tiepoint_accuracy_pix = args["reference_args"]["tiepoint_accuracy_pix"]
        chunk.camera_location_accuracy = Metashape.Vector(
            [camera_accuracy_m, camera_accuracy_m, camera_accuracy_m])
        chunk.camera_rotation_accuracy = Metashape.Vector(
            [camera_accuracy_deg, camera_accuracy_deg, camera_accuracy_deg])
        chunk.marker_location_accuracy = Metashape.Vector(
            [marker_accuracy_m, marker_accuracy_m, marker_accuracy_m])
        chunk.marker_projection_accuracy = float(marker_accuracy_pix)
        chunk.tiepoint_accuracy = float(tiepoint_accuracy_pix)

        document.save()
    return document


def accessMaxChunk(document, chunk_name):
    # [2]  access chunk
    if chunk_name:
        # choose the chunk with given name
        for chunk_i in document.chunks:
            if not chunk_i.enabled:
                continue
            if chunk_i.label == chunk_name:
                return chunk_i
    else:
        # choose the chunk with maximum images
        chunk_images_num = 0
        for chunk_i in document.chunks:
            if not chunk_i.enabled:
                continue
            chunki_images_num = len(chunk_i.cameras)
            if chunki_images_num >= chunk_images_num:
                chunk_images_num = chunki_images_num
                chunk = chunk_i
            else:
                continue
        return chunk


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
    if metashape_version == 1.8:
        chunk.optimizeCameras(fit_f=args['fit_f'],
                              fit_cx=args['fit_cx'], fit_cy=args['fit_cy'],
                              fit_b1=args['fit_b1'], fit_b2=args['fit_b2'],
                              fit_k1=args['fit_k1'], fit_k2=args['fit_k2'],
                              fit_k3=args['fit_k3'], fit_k4=args['fit_k4'],
                              fit_p1=args['fit_p1'], fit_p2=args['fit_p2'],
                              fit_corrections=args['fit_corrections'],
                              adaptive_fitting=args['adaptive_fitting'],
                              tiepoint_covariance=args['tiepoint_covariance'])
    elif metashape_version == 1.5:
        chunk.optimizeCameras(fit_f=args['fit_f'],
                              fit_cx=args['fit_cx'], fit_cy=args['fit_cy'],
                              fit_b1=args['fit_b1'], fit_b2=args['fit_b2'],
                              fit_k1=args['fit_k1'], fit_k2=args['fit_k2'],
                              fit_k3=args['fit_k3'], fit_k4=args['fit_k4'],
                              fit_p1=args['fit_p1'], fit_p2=args['fit_p2'],
                              fit_corrections=args['fit_corrections'],
                              adaptive_fitting=args['adaptive_fitting'],
                              tiepoint_covariance=args['tiepoint_covariance'])


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
    if metashape_version == 1.8:
        chunk.buildOrthomosaic(surface_data=Metashape.ElevationData, blending_mode=Metashape.MosaicBlending,
                               fill_holes=True,
                               ghosting_filter=True, cull_faces=False, refine_seamlines=False)
    elif metashape_version == 1.5:
        chunk.buildOrthomosaic(surface=Metashape.ElevationData, blending=Metashape.MosaicBlending, fill_holes=True,
                               cull_faces=False, refine_seamlines=False)


def buildTiledModel(chunk, TiledModel_pixel_size, TiledModel_face_count):
    if metashape_version == 1.8:
        chunk.buildTiledModel(source_data=Metashape.DataSource.DenseCloudData,
                              pixel_size=TiledModel_pixel_size, tile_size=256, face_count=TiledModel_face_count,
                              ghosting_filter=False, transfer_texture=False, keep_depth=True, merge=False)
    elif metashape_version == 1.5:
        chunk.buildTiledModel(tile_size=256, face_count=TiledModel_face_count,
                              ghosting_filter=False, transfer_texture=False, keep_depth=True)


def exportDSM(chunk, path_to_save):
    path_DSM = path_to_save + "_DSM.tif"
    if metashape_version == 1.8:
        chunk.exportRaster(path=path_DSM, source_data=Metashape.ElevationData, nodata_value=-32767,
                           save_kml=False, save_world=True, save_scheme=False, save_alpha=False,
                           image_compression=image_compression_paras,
                           image_description=os.path.basename(path_to_save), white_background=False,
                           title='Digital Surface Model')


def exportDOM(chunk, path_to_save):
    path_DOM = path_to_save + "_DOM.tif"
    if metashape_version == 1.8:
        chunk.exportRaster(path=path_DOM, source_data=Metashape.OrthomosaicData, nodata_value=-32767,
                           save_kml=False, save_world=True, save_scheme=False, save_alpha=False,
                           image_compression=image_compression_paras,
                           image_description=os.path.basename(path_to_save), white_background=False,
                           title='Orthomosaic')


def exportTiledModel(chunk, path_to_save, TiledModel_TiledModelFormat, TiledModel_ModelFormat, TiledModel_ImageFormat):
    path_TM = path_to_save + f"_TiledModel_{TiledModel_TiledModelFormat}.zip"
    if metashape_version == 1.8:
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
