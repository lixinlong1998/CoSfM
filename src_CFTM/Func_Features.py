from config import *
import cv2
import json
import numpy as np
import src_CFTM.Func_Files as Wrif
import src_CFTM.ConnectData_COLMAP as Colmap


# print(f'OpenCV: {cv2.__version__} for python installed and working')  # 4.7.0


def extractFeature(args):
    feature_extraction_mode = args["feature_extraction_mode"]
    data_package_path = args["data_package_path"]
    feature_path = args["feature_path"]
    Cameras = Wrif.readDictionary_Int_List(data_package_path + '/Cameras.txt')
    cameraPaths = Wrif.readDictionary_Str_Int(data_package_path + '/cameraPaths.txt')

    if feature_extraction_mode == 0:
        # 0: use OpenCV to extract SIFT feature without CUDA (default but slow);
        chunkKeypoints, chunkDescriptors = extractSIFT(Cameras)
    elif feature_extraction_mode == 1:
        # 1: use OpenCV to extract SIFT-GPU feature with CUDA (need additional manual compilation);
        chunkKeypoints, chunkDescriptors = extractSIFT_CUDA(Cameras, cameraPaths, feature_path)
    elif feature_extraction_mode == 2:
        # 2: use COLMAP to extract SIFT-GPU feature (recommended).
        chunkKeypoints, chunkDescriptors = Colmap.extractSIFT(Cameras, cameraPaths, feature_path)

    # save the keypoints and descriptors
    # with open(os.path.join(data_package_path, 'chunkKeypoints.json'), 'w') as json_file:
    #     json.dump(chunkKeypoints, json_file)
    # with open(os.path.join(data_package_path, 'chunkDescriptors.json'), 'w') as json_file:
    #     json.dump(chunkDescriptors, json_file)
    np.savez(os.path.join(data_package_path, 'chunkKeypoints.npz'), **chunkKeypoints)
    np.savez(os.path.join(data_package_path, 'chunkDescriptors.npz'), **chunkDescriptors)
    return None


def extractSIFT(Cameras):
    '''
    input:
        Cameras = {camera_id: Camera, ...}
            Camera = [transform, reference, covariance, sensors_id, state, path, identifier]
                transform = [camera_transform, camera_loc, camera_rot, camera_transform_PJCS]
                reference = [camera_loc_Ref, location_accuracy, camera_rot_Ref, rotation_accuracy]
                covariance = [camera_location_covariance, camera_rotation_covariance]
                states = [label, camera_loc_Ref, camera_rot_Ref, selected, camera_orientation]
    output:
        chunkKeypoints = {identifier: 2darray.shape(numKP, 6).dtype(float32), ...}
        chunkDescriptors = {identifier: 2darray.shape(numKP, 128).dtype(uint8), ...}
    '''
    chunkKeypoints = {}
    chunkDescriptors = {}
    for camera_id, Camera in Cameras.items():
        identifier = Camera[-1]
        image_path = Camera[5]
        # Read the image
        image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
        # Initialize the SIFT detector
        sift = cv2.SIFT_create()
        # Detect keypoints and compute descriptors
        keypoints, descriptors = sift.detectAndCompute(image, None)
        # Store information about keypoints in an array
        keypoints_info = np.array([(kp.pt[0], kp.pt[1], kp.size, kp.angle, kp.response) for kp in keypoints],
                                  dtype=np.float32)
        # Store keypoints and descriptors in a dictionary
        chunkKeypoints[str(identifier)] = keypoints_info
        chunkDescriptors[str(identifier)] = descriptors
    return dict(chunkKeypoints), dict(chunkDescriptors)


def extractSIFT_CUDA(Cameras, cameraPaths, feature_path):
    '''
    input:
        Cameras = {camera_id: Camera, ...}
            Camera = [transform, reference, covariance, sensors_id, state, path, identifier]
                transform = [camera_transform, camera_loc, camera_rot, camera_transform_PJCS]
                reference = [camera_loc_Ref, location_accuracy, camera_rot_Ref, rotation_accuracy]
                covariance = [camera_location_covariance, camera_rotation_covariance]
                states = [label, camera_loc_Ref, camera_rot_Ref, selected, camera_orientation]
    output:
        chunkKeypoints = {identifier: 2darray.shape(numKP, 6).dtype(float32), ...}
        chunkDescriptors = {identifier: 2darray.shape(numKP, 128).dtype(uint8), ...}
    '''
    chunkKeypoints = {}
    chunkDescriptors = {}
    for camera_id, Camera in Cameras.items():
        identifier = Camera[-1]
        image_path = Camera[5]
        # Read the image
        image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
        # Create a SIFT object
        sift = cv2.cuda.SIFT_create()
        # Upload the image to the GPU
        d_image = cv2.cuda_GpuMat()
        d_image.upload(image)
        # Perform SIFT computation on the GPU
        d_keypoints, d_descriptors = sift.detectAndComputeAsync(d_image, None)
        # Download keypoints and descriptors from the GPU to the host
        keypoints = d_keypoints.download()
        descriptors = d_descriptors.download()
        # Store information about keypoints in an array
        keypoints_info = np.array([(kp.pt[0], kp.pt[1], kp.size, kp.angle, kp.response) for kp in keypoints],
                                  dtype=np.float32)
        # Store keypoints and descriptors in a dictionary
        chunkKeypoints[str(identifier)] = keypoints_info
        chunkDescriptors[str(identifier)] = descriptors
    return dict(chunkKeypoints), dict(chunkDescriptors)
