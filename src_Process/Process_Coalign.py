from config import *
import Metashape
from src_CFTM import Func_Files
from src_CFTM import Func_CommonTiePoints
from src_Metashape import FuncMs_Basic


def splitChunk(doc, coalign_chunk_name, epochs, epoch_mode):
    '''
    function:
        split the coalign chunk into individual chunks with its own epoch images
    parameters should be given as follows:
        coalign_chunk
        epochs = [['2020:10:14'],['2021:07:03'],['2021:12:03']]
    return:
        chunk_1,chunk_2,chunk_3
    '''
    coalign_chunk = FuncMs_Basic.accessMaxChunk(doc, coalign_chunk_name)

    chunksNum = len(epochs)
    chunk_name_list = [f'chunk_{i}' for i in range(chunksNum)]

    # check if there already have splited chunks
    chunk_already_exist = []
    for chunk_i in doc.chunks:
        if chunk_i.label in chunk_name_list:
            chunk_already_exist.append(chunk_i)
    if len(chunk_already_exist) == len(chunk_name_list):
        # there already have splited chunks
        return chunk_already_exist

    # split coalign chunk into chunks
    else:
        chunk_list = []
        for chunk_name in chunk_name_list:
            locals()[chunk_name] = coalign_chunk.copy(keypoints=False)
            locals()[chunk_name].label = chunk_name
            chunk_list.append(locals()[chunk_name])
        # remove cameras which capture date is not in line with the epoch of chunk
        for i, chunk in enumerate(chunk_list):
            # remove cameras which capture date is not in line with the epoch of chunk
            camerasRemoved_list = []
            for camera in chunk.cameras:
                # check which epoch dose this photo belongs to
                if epoch_mode == "DATE":
                    # get the capture date of photo
                    Date = camera.photo.meta['Exif/DateTimeOriginal'].split(' ')[0]  # return '2021:04:04'
                    if Date not in epochs[i]:
                        camerasRemoved_list.append(camera)
                elif epoch_mode == "FOLDER":
                    camera_path = camera.photo.path
                    if camera_path != epochs[i]:
                        camerasRemoved_list.append(camera)
            chunk.remove(camerasRemoved_list)
        return chunk_list


def reconstruction(doc, args):
    starttime0 = time.perf_counter()

    # [0]  open document
    process_reconstruction = args["process_args"]["process_reconstruction"]
    process_build_DSM = args["process_args"]["process_build_DSM"]
    process_build_DOM = args["process_args"]["process_build_DOM"]
    process_build_TiledModel = args["process_args"]["process_build_TiledModel"]
    process_export_DSM = args["process_args"]["process_export_DSM"]
    process_export_DOM = args["process_args"]["process_export_DOM"]
    process_export_TiledModel = args["process_args"]["process_export_TiledModel"]
    process_white_balance = args["process_args"]["process_white_balance"]
    save_project_each_step = args["process_args"]["save_project_each_step"]
    quality_DPC = args["photogrammetry_args"]["quality_DPC"]
    dense_cloud_FileterMode = args["photogrammetry_args"]["dense_cloud_FileterMode"]
    tiled_model_TiledModelFormat = args["photogrammetry_args"]["tiled_model_TiledModelFormat"]
    tiled_model_pixel_size = args["photogrammetry_args"]["tiled_model_pixel_size"]
    tiled_model_face_count = args["photogrammetry_args"]["tiled_model_face_count"]
    tiled_model_ModelFormat = args["photogrammetry_args"]["tiled_model_ModelFormat"]
    tiled_model_ImageFormat = args["photogrammetry_args"]["tiled_model_ImageFormat"]
    data_package_path = args["data_package_path"]
    workspace_path = args["workspace_path"]
    project_name = args["project_name"]
    coalign_chunk_name = args["coalign_chunk_name"]
    epoch_mode = args["epoch_mode"]
    camera_ids = Func_Files.readDictionary_Int_Int(data_package_path + '/camera_ids.txt')
    cameraPaths = Func_Files.readDictionary_Str_Int(data_package_path + '/cameraPaths.txt')
    if epoch_mode == "DATE":
        epochs = Func_CommonTiePoints.getEpoch_byDATE(camera_ids)
    elif epoch_mode == "FOLDER":
        epochs = Func_CommonTiePoints.getEpoch_byFOLDER(cameraPaths)

    # [1]  split coaligned chunk into individual chunks
    chunk_list = splitChunk(doc, coalign_chunk_name, epochs, epoch_mode)
    doc.save()

    for chunk_id, chunk in enumerate(chunk_list):
        print(chunk.label)
        path_to_save = workspace_path + f"/{os.path.splitext(project_name)[0]}_{chunk_id + 1}"

        # [2]  build dense cloud for each individual chunks
        if process_reconstruction and not chunk.dense_cloud:
            FuncMs_Basic.buildDenseCloud(chunk, quality_DPC, dense_cloud_FileterMode)
            if save_project_each_step:
                doc.save()

        # [3]  build DSM
        if process_build_DSM and not chunk.elevation:
            FuncMs_Basic.buildDSM(chunk)
            if save_project_each_step:
                doc.save()

        # [4]  calibrate colors and white balance
        if process_white_balance and chunk.elevation and not chunk.orthomosaic:
            chunk.calibrateColors(source_data=Metashape.DataSource.ElevationData, white_balance=True)
            if save_project_each_step:
                doc.save()

        # [5] build DOM
        if process_build_DOM and chunk.elevation and not chunk.orthomosaic:
            FuncMs_Basic.buildDOM(chunk)
            if save_project_each_step:
                doc.save()

        # [6] build TiledModel
        if process_build_TiledModel and chunk.elevation and not chunk.tiled_model:
            FuncMs_Basic.buildTiledModel(chunk, tiled_model_pixel_size, tiled_model_face_count)
            if save_project_each_step:
                doc.save()

        # [7] export DSM, DOM and TiledModel
        if chunk.elevation and process_export_DSM:
            FuncMs_Basic.exportDSM(chunk, path_to_save)
        if chunk.orthomosaic and process_export_DOM:
            FuncMs_Basic.exportDOM(chunk, path_to_save)
        if chunk.tiled_model and process_export_TiledModel:
            FuncMs_Basic.exportTiledModel(chunk, path_to_save, tiled_model_TiledModelFormat, tiled_model_ModelFormat,
                                          tiled_model_ImageFormat)

        # [8]  export report
        FuncMs_Basic.exportReport(chunk, path_to_save)

        # save project
        doc.save()
    print('[Script]    [Script][TimeCost]    Reconstruction complete in :', time.perf_counter() - starttime0)
