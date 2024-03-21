from config import *
import subprocess
import sqlite3
import numpy as np

MAX_IMAGE_ID = 2147483647  # 2 ** 31 - 1


def extractSIFT(Cameras, cameraPaths, feature_path):
    '''
    input:
        Cameras = {camera_id: Camera, ...}
            Camera = [transform, reference, covariance, sensors_id, state, path, identifier]
                transform = [camera_transform, camera_loc, camera_rot, camera_transform_PJCS]
                reference = [camera_loc_Ref, location_accuracy, camera_rot_Ref, rotation_accuracy]
                covariance = [camera_location_covariance, camera_rotation_covariance]
                states = [label, camera_loc_Ref, camera_rot_Ref, selected, camera_orientation]
        cameraPaths = {camera.photo.path: identifier, ...}
    output:
        chunkKeypoints = {identifier: 2darray.shape(numKP, 6).dtype(float32), ...}
        chunkDescriptors = {identifier: 2darray.shape(numKP, 128).dtype(uint8), ...}
    '''
    # Collect all parent folder paths of the images
    colmap_workspace = creat_colmap_workspace(Cameras, feature_path)

    # Create a COLMAP project file for each image folder and perform feature extraction
    for colmap_path_dir in colmap_workspace:
        # load paths
        image_folder_path = colmap_path_dir["image_folder_path"]
        colmap_database_path = colmap_path_dir["colmap_database_path"]
        colmap_project_path = colmap_path_dir["colmap_project_path"]
        colmap_output_path = colmap_path_dir["colmap_output_path"]

        # set project argument
        colmap_quality = 'high'  # {low, medium, high, extreme}

        print('project:', colmap_project_path)
        print('database:', colmap_database_path)
        print('image folder:', image_folder_path)

        # # create project
        # project_generator(colmap_project_path, colmap_output_path, colmap_quality)

        # # create database
        # database_creator(colmap_project_path, colmap_database_path)

        # extract feature in colmap
        feature_extractor(colmap_project_path, colmap_database_path, image_folder_path)

    # Collect feature results
    chunkKeypoints, chunkDescriptors = getFeatures(colmap_workspace, cameraPaths)
    return chunkKeypoints, chunkDescriptors


def extractMatch(Cameras, cameraPaths, feature_path):
    '''
    input:
        Cameras = {camera_id: Camera, ...}
            Camera = [transform, reference, covariance, sensors_id, state, path, identifier]
                transform = [camera_transform, camera_loc, camera_rot, camera_transform_PJCS]
                reference = [camera_loc_Ref, location_accuracy, camera_rot_Ref, rotation_accuracy]
                covariance = [camera_location_covariance, camera_rotation_covariance]
                states = [label, camera_loc_Ref, camera_rot_Ref, selected, camera_orientation]
        cameraPaths = {camera.photo.path: identifier, ...}
    output:
        chunkKeypoints = {identifier: 2darray.shape(numKP, 6).dtype(float32), ...}
        chunkDescriptors = {identifier: 2darray.shape(numKP, 128).dtype(uint8), ...}
    '''

    # Collect all parent folder paths of the images
    colmap_workspace = creat_colmap_workspace(Cameras, feature_path)

    # Create a COLMAP project file for each image folder and perform feature extraction
    for colmap_path_dir in colmap_workspace:

        # load paths
        colmap_database_path = colmap_path_dir["colmap_database_path"]
        colmap_project_path = colmap_path_dir["colmap_project_path"]

        if os.path.exists(colmap_database_path):
            exhaustive_matcher(colmap_project_path, colmap_database_path)
        else:
            raise Exception(f'Database: {colmap_database_path} is not found!')

    # Collect match results
    chunkMatches = getMatches(colmap_workspace, cameraPaths, min_num_matches=0)
    return chunkMatches


def extractInlierMatch(Cameras, cameraPaths, feature_path):
    '''
    input:
        Cameras = {camera_id: Camera, ...}
            Camera = [transform, reference, covariance, sensors_id, state, path, identifier]
                transform = [camera_transform, camera_loc, camera_rot, camera_transform_PJCS]
                reference = [camera_loc_Ref, location_accuracy, camera_rot_Ref, rotation_accuracy]
                covariance = [camera_location_covariance, camera_rotation_covariance]
                states = [label, camera_loc_Ref, camera_rot_Ref, selected, camera_orientation]
        cameraPaths = {camera.photo.path: identifier, ...}
    output:
        chunkKeypoints = {identifier: 2darray.shape(numKP, 6).dtype(float32), ...}
        chunkDescriptors = {identifier: 2darray.shape(numKP, 128).dtype(uint8), ...}
    '''

    # Collect all parent folder paths of the images
    colmap_workspace = creat_colmap_workspace(Cameras, feature_path)

    # Collect match results
    chunkInlierMatches, chunkTwoViewsGeometry = getInlierMatches(colmap_workspace, cameraPaths,
                                                                 min_num_matches=0, geometry=True)
    return chunkInlierMatches, chunkTwoViewsGeometry


########################################################################################################################
def project_generator(colmap_project_path, colmap_output_path, colmap_quality):
    colmap_command = [
        PATH_COLMAP_BAT, 'project_generator',
        '--project_path', colmap_project_path,
        '--output_path', colmap_output_path,
        '--quality', colmap_quality,
    ]
    subprocess.run(colmap_command)


def database_creator(colmap_project_path, colmap_database_path):
    colmap_command = [
        PATH_COLMAP_BAT, 'project_generator',
        '--project_path', colmap_project_path,
        '--database_path', colmap_database_path,
    ]
    subprocess.run(colmap_command)


def feature_extractor(colmap_project_path, colmap_database_path, image_folder_path):
    colmap_command = [
        PATH_COLMAP_BAT, 'feature_extractor',
        # '--random_seed ', '0',
        # '--log_to_stderr ', '0',
        # '--log_level ', '2',
        #
        # '--project_path', colmap_project_path,
        '--database_path', rf"{colmap_database_path}",
        '--image_path', rf"{image_folder_path}",
        # '--image_list_path arg
        # '--ImageReader.mask_path arg

        # '--ImageReader.camera_model ', 'SIMPLE_RADIAL',
        # '--ImageReader.single_camera ', '0',
        # '--ImageReader.single_camera_per_folder ', '0',
        # '--ImageReader.existing_camera_id ', '-1',
        # '--ImageReader.camera_params arg
        # '--ImageReader.default_focal_length_factor ', '1.2',
        # '--ImageReader.camera_mask_path arg
        #
        # '--SiftExtraction.num_threads ', '-1',
        # '--SiftExtraction.use_gpu ', '1',
        # '--SiftExtraction.gpu_index ', '-1',
        # '--SiftExtraction.max_image_size ', '3200',
        # '--SiftExtraction.max_num_features ', '8192',
        # '--SiftExtraction.first_octave ', '-1',
        # '--SiftExtraction.num_octaves ', '4',
        # '--SiftExtraction.octave_resolution ', '3',
        # '--SiftExtraction.peak_threshold ', '0.0066666666666666671',
        # '--SiftExtraction.edge_threshold ', '10',
        # '--SiftExtraction.estimate_affine_shape ', '0',
        # '--SiftExtraction.max_num_orientations ', '2',
        # '--SiftExtraction.upright ', '0',
        # '--SiftExtraction.domain_size_pooling ', '0',
        # '--SiftExtraction.dsp_min_scale ', '0.16666666666666666',
        # '--SiftExtraction.dsp_max_scale ', '3',
        # '--SiftExtraction.dsp_num_scales ', '10'
    ]
    # print(colmap_command)
    subprocess.run(colmap_command)


def exhaustive_matcher(colmap_project_path, colmap_database_path):
    colmap_command = [
        PATH_COLMAP_BAT, 'exhaustive_matcher',
        # '--random_seed', '0',
        # '--log_to_stderr', '0',
        # '--log_level', '2',

        '--project_path', os.path.normpath(colmap_project_path),
        '--database_path', os.path.normpath(colmap_database_path),

        '--SiftMatching.num_threads', '-1',
        '--SiftMatching.use_gpu', '1',
        '--SiftMatching.gpu_index', '-1',
        '--SiftMatching.max_ratio', '0.80000000000000004',
        '--SiftMatching.max_distance', '0.69999999999999996',
        '--SiftMatching.cross_check', '1',
        '--SiftMatching.max_error', '4',
        '--SiftMatching.max_num_matches', '32768',
        '--SiftMatching.confidence', '0.999',
        '--SiftMatching.max_num_trials', '10000',
        '--SiftMatching.min_inlier_ratio', '0.25',
        '--SiftMatching.min_num_inliers', '15',
        '--SiftMatching.multiple_models', '0',
        '--SiftMatching.guided_matching', '0',
        '--ExhaustiveMatching.block_size', '50'
    ]
    subprocess.run(colmap_command)


def getFeatures(colmap_workspace, cameraPaths):
    '''
    input:
        colmap_workspace = [colmap_path_dir,...]
        colmap_path_dir = {
            "image_folder_path": image_folder_path,
            "colmap_folder_path": colmap_folder_path,
            "colmap_database_path": os.path.join(colmap_folder_path, f'database_{image_folder_id}.db'),
            "colmap_project_path": os.path.join(colmap_folder_path, f'project_{image_folder_id}.ini'),
            "colmap_output_path": os.path.join(colmap_folder_path, f'output_{image_folder_id}'),
        }
        cameraPaths = {camera.photo.path : identifier,...}
    output:
        chunkKeypoints = {identifier:2darray.shape(numKP,6).dtype(float32),...}
            ,where 2darray = [[u,v,affinity1, affinity2, affinity3, affinity4],...]
        chunkDescriptors = {identifier:2darray.shape(numKP,128).dtype(uint8),...}
    '''
    chunkKeypoints = []
    chunkDescriptors = []
    for colmap_path_dir in colmap_workspace:

        # load paths
        image_folder_path = colmap_path_dir["image_folder_path"]
        colmap_database_path = colmap_path_dir["colmap_database_path"]
        print(image_folder_path)
        print(colmap_database_path)

        # Load the database
        database = sqlite3.connect(colmap_database_path)
        Images = dict((image_id, [name, camera_id])
                      for image_id, name, camera_id in
                      database.execute("SELECT image_id, name, camera_id FROM images"))
        Keypoints = dict((image_id, blob_to_array(data, np.float32, (-1, 6)))
                         for image_id, data in database.execute("SELECT image_id, data FROM keypoints"))
        Descriptors = dict((image_id, blob_to_array(data, np.uint8, (-1, 128)))
                           for image_id, data in database.execute("SELECT image_id, data FROM descriptors"))

        # Replace the index with the identifier
        image_folder_keypoints = {}
        image_folder_descriptors = {}
        for image_id, value in Images.items():
            image_name = value[0]
            image_path = image_folder_path + '/' + image_name
            identifier = cameraPaths[image_path]
            image_folder_keypoints[str(identifier)] = Keypoints[image_id]  # np.savez() keywords must be strings
            image_folder_descriptors[str(identifier)] = Descriptors[image_id]
        chunkKeypoints += list(image_folder_keypoints.items())
        chunkDescriptors += list(image_folder_descriptors.items())
        database.close()
    return dict(chunkKeypoints), dict(chunkDescriptors)


def getMatches(colmap_workspace, cameraPaths, min_num_matches=0):
    '''
    input:
        colmap_workspace = [colmap_path_dir,...]
        colmap_path_dir = {
            "image_folder_path": image_folder_path,
            "colmap_folder_path": colmap_folder_path,
            "colmap_database_path": os.path.join(colmap_folder_path, f'database_{image_folder_id}.db'),
            "colmap_project_path": os.path.join(colmap_folder_path, f'project_{image_folder_id}.ini'),
            "colmap_output_path": os.path.join(colmap_folder_path, f'output_{image_folder_id}'),
        }
        cameraPaths = {camera.photo.path : identifier,...}
    output:
        chunkMatches = {(identifier,identifier):2darray.shape(numMatch,2).dtype(np.uint32),...}
            ,where 2darray = [[keypoint_id1,keypoint_id2],...]
    '''
    chunkMatches = []
    for colmap_path_dir in colmap_workspace:

        # load paths
        image_folder_path = colmap_path_dir["image_folder_path"]
        colmap_database_path = colmap_path_dir["colmap_database_path"]

        # load database
        database = sqlite3.connect(colmap_database_path)
        cursor = database.cursor()

        # get images
        Images = dict((image_id, [name, camera_id])
                      for image_id, name, camera_id in
                      database.execute("SELECT image_id, name, camera_id FROM images"))

        # get matches
        Matches = {}
        cursor.execute("SELECT pair_id, data, rows FROM Matches WHERE rows>=?;", (min_num_matches,))
        for row in cursor:

            pair_id = row[0]
            # connect index to identifier
            image_id1, image_id2 = pair_id_to_image_ids(pair_id)
            image_name1, image_name2 = Images[image_id1][0], Images[image_id2][0]
            identifier1 = cameraPaths[os.path.join(image_folder_path, image_name1)]
            identifier2 = cameraPaths[os.path.join(image_folder_path, image_name2)]

            if not row[1] or row[2] == 0:
                continue
            else:
                inlier_matches = np.frombuffer(row[1], dtype=np.uint32).reshape(-1, 2)
                Matches[(identifier1, identifier2)] = inlier_matches
        chunkMatches += list(Matches.items())

        cursor.close()
        database.close()
    return dict(chunkMatches)


def getInlierMatches(colmap_workspace, cameraPaths, min_num_matches=0, geometry=True):
    '''
    input:
        colmap_workspace = [colmap_path_dir,...]
        colmap_path_dir = {
            "image_folder_path": image_folder_path,
            "colmap_folder_path": colmap_folder_path,
            "colmap_database_path": os.path.join(colmap_folder_path, f'database_{image_folder_id}.db'),
            "colmap_project_path": os.path.join(colmap_folder_path, f'project_{image_folder_id}.ini'),
            "colmap_output_path": os.path.join(colmap_folder_path, f'output_{image_folder_id}'),
        }
        cameraPaths = {camera.photo.path : identifier,...}
    output:
        chunkInlierMatches = {(identifier,identifier):2darray.shape(numMatch,2).dtype(np.uint32),...}
        chunkTwoViewsGeometry = {(identifier,identifier):[F,E,H]],...}
            ,where F,E,H = 2darray.shape(3,3).dtype(np.float64)
    '''
    chunkInlierMatches = []
    chunkTwoViewsGeometry = []
    for colmap_path_dir in colmap_workspace:

        # load paths
        image_folder_path = colmap_path_dir["image_folder_path"]
        colmap_database_path = colmap_path_dir["colmap_database_path"]

        # load database
        database = sqlite3.connect(colmap_database_path)
        cursor = database.cursor()

        # get images
        Images = dict((image_id, [name, camera_id])
                      for image_id, name, camera_id in
                      database.execute("SELECT image_id, name, camera_id FROM images"))
        # get InlierMatches
        InlierMatches = {}
        TwoViewGeometries = {}
        cursor.execute("SELECT pair_id, data,rows,F,E,H FROM two_view_geometries WHERE rows>=?;", (min_num_matches,))
        for row in cursor:
            pair_id = row[0]
            # connect index to identifier
            image_id1, image_id2 = pair_id_to_image_ids(pair_id)
            image_name1, image_name2 = Images[image_id1][0], Images[image_id2][0]
            identifier1 = cameraPaths[os.path.join(image_folder_path, image_name1)]
            identifier2 = cameraPaths[os.path.join(image_folder_path, image_name2)]
            # get InlierMatches and TwoViewGeometries
            if not row[1] or row[2] == 0:
                continue
            else:
                inlier_matches = np.frombuffer(row[1], dtype=np.uint32).reshape(-1, 2)
                InlierMatches[(identifier1, identifier2)] = inlier_matches
                if geometry:
                    F = np.fromstring(row[3], dtype=np.float64).reshape(3, 3)
                    E = np.fromstring(row[4], dtype=np.float64).reshape(3, 3)
                    H = np.fromstring(row[5], dtype=np.float64).reshape(3, 3)
                    TwoViewGeometries[(identifier1, identifier2)] = [F, E, H]
        chunkInlierMatches += list(InlierMatches.items())
        if geometry:
            chunkTwoViewsGeometry += list(TwoViewGeometries.items())
        cursor.close()
        database.close()
    return dict(chunkInlierMatches), dict(chunkTwoViewsGeometry)


########################################################################################################################

def creat_colmap_workspace(Cameras, feature_path):
    # Collect all parent folder paths of the images
    image_path_list = [Camera[5] for Camera in Cameras.values()]
    image_folder_set = set()
    for image_path in image_path_list:
        image_folder_set.add(os.path.dirname(image_path))
    image_folder_list = list(image_folder_set)

    # Create a COLMAP project file for each image folder
    colmap_workspace = []
    for image_folder_path in image_folder_list:
        image_folder_name = os.path.basename(image_folder_path)
        colmap_folder_path = os.path.join(feature_path, f'colmap_{image_folder_name}')
        os.makedirs(colmap_folder_path, exist_ok=True)

        colmap_path_dir = {
            "image_folder_path": image_folder_path,
            "colmap_folder_path": colmap_folder_path,
            "colmap_database_path": os.path.join(colmap_folder_path, f'database_{image_folder_name}.db'),
            "colmap_project_path": os.path.join(colmap_folder_path, f'project_{image_folder_name}.ini'),
            "colmap_output_path": os.path.join(colmap_folder_path, f'output_{image_folder_name}'),
        }

        os.makedirs(colmap_path_dir["colmap_output_path"], exist_ok=True)
        colmap_workspace.append(colmap_path_dir)

    return colmap_workspace


def blob_to_array(blob, dtype, shape=(-1,)):
    return np.fromstring(blob, dtype=dtype).reshape(*shape)


def array_to_blob(array):
    return array.tostring()


def image_ids_to_pair_id(image_id1, image_id2):
    if image_id1 > image_id2:
        image_id1, image_id2 = image_id2, image_id1
    return image_id1 * MAX_IMAGE_ID + image_id2


def pair_id_to_image_ids(pair_id):
    image_id2 = pair_id % MAX_IMAGE_ID
    # image_id1 = (pair_id - image_id2) / MAX_IMAGE_ID
    image_id1 = pair_id // MAX_IMAGE_ID
    return image_id1, image_id2
