from config import *
import Metashape
from src_CFTM import Func_Breakpoint, Func_Features, Func_Files
from src_Metashape import FuncMs_Basic
from src_Process import Process_Analysis
from src_Process import Process_Coalign
from src_Process import Process_CTPsGeneration_part1
from src_Process import Process_IterativeOptimization_part1
from src_Process import Process_IterativeOptimization_part2
from src_Process import Process_IterativeOptimization_part3

'''CoSfM version 1.0

source: https://github.com/LXLone/CoSfM
author: Xinlong Li

# Introduction
The Common Tie Points(CTPs) generator consists of three parts:
    Part1: Export and store data from Metashape;
    Part2: Read the data and generate CTPs by using external functions in 'src_CFTM', then save the results as .txt files;
    Part3: Read the files and import CTPs with Metashape.Marker format into Metashape, then run the optimization;

# Preparation
Highly recommend checking the README or the user manual book before use!
Here is a short check list for setting up:
    1. Add the Metashape into your system environment variable of PATH;
    2. Install the extra third libraries;
    3. Set up the paths in ./config.py;
    4. Set up the parameters in SETUP section below in this script;
    5. The docment has been created and the images have been merged into the same chunk with proper setting.

# Run script
Since the script contains parallel computing, you need to run this script from cmd by using the following command:

metashape.exe -r G:\AResearchG\20221223_CoSfM\Release\CFTM_v1.0\Run.py

(please change the root of script to yours.)
'''
########################################################################################################################
################################################       SETUP       #####################################################
# Open a project file
workspace_path = r"K:\CoSfM_v1\project_with_cftm"
project_name = "cftm_test_project_test2.psx"
images_path = r'K:\CoSfM_v1\images'
check_points_path = r"K:\CoSfM_v1/cftm_test_project_CPsdatabase.txt"

# please make sure your major version is:
compatible_major_version = "1.8"

# Process arguments
process_args = {
    # Choose the processing options for CFTM.
    "process_coalign": False,
    "process_feature_extraction": True,
    "process_CTPs_generation": True,
    "process_iterative_optimization": True,
    "process_analysis": True,
    "process_analysis_matches": False,
    "process_compare_matches": False,

    # Choose the processing options for reconstruction.
    "process_reconstruction": True,
    "process_build_DSM": True,
    "process_build_DOM": False,
    "process_build_TiledModel": False,
    "process_export_DSM": True,
    "process_export_DOM": False,
    "process_export_TiledModel": False,
    "process_white_balance": False,

    "save_project_each_step": True,
}

# Select the mode for dividing epochs.
# If it's "DATE," It will divide them by day based on the dates of the photos.
# If it's "FOLDER," it will divide them based on the folders where the photos are stored.
epoch_mode = 'DATE'  # or 'FOLDER'

# Which feature extraction mode do you want to use?
# 0 : use OpenCV to extract SIFT feature without cuda (default but slow);
# 1 : use OpenCV to extract SIFT-GPU feature with cuda (need additional manual compilation);
# 2 : use COLMAP to extract SIFT-GPU feature (recommend).
feature_extraction_mode = 2

# Whether to continue the run from the last breakpoint?
begin_from_breakpoint = True

# Specify which chunk to process; leave empty will automatically select the chunk which has maximum images.
coalign_chunk_name = ''

# Generate CTPs for all grids or only for those grid without CTPs.
for_all_grids = True

# If CTPs are generated for all grids, do you what to reuse the CTPs generated from co-alignment?
reuse_ICTPs = True

# Whether to use *.shp to mask the unstable areas?
unstable_areas_mask_path = ''  # 'I:/20230418_CTPsGenerator/Adata/Region/UnstableRegion.shp'

# Common Feature Track Matching
'''helper:
gird_size_1:        	Divide the point cloud of the survey area into squares with given size (in meters) for CFTM.
gird_size_2:        	Divide the point cloud of the survey area into squares with given size (in meters) for iterative optimization.
scoring:            	This is a spatial uniformity scoring table, with a pyramid grid constructed from the outer rectangle of the survey area,
                    	where the length of the first square is 1 metre, the length of the second square is 2 metres, and the length of the nth square is 2^(n-1) metres;
                    	whenever a point drops into a square of a certain layer, the total score is added based on the scoring table.
pool_size:          	Number of CPU cores enabled.
criterias_1:        	Maximum reprojection error threshold and minimum intersection Angle threshold in feature matching stage, in pixels and degree respectively.
criterias_2:        	Maximum reprojection error threshold and minimum intersection Angle threshold in tracks filtering stage, in pixels and degree respectively.
ratio_inlier_init:  	An initial guess of the average ratio of inliers in track, in the range from 0 to 1.
confidence:         	A probability of being able to get the consensus set by randomly sample certain times, in the range from 0 to 1.
threshold_RepError: 	A threshold used for the final assessment of point quality of track, in pixels.
threshold_Distance: 	For track matching, only track pairs whose distance between points is less than a threshold are considered, in meters.
num_nominated:      	Number of pre-selected CTPs based on their quality.
num_selected:    		Number of the final choosed CTPs based on their quality.
num_max_iterations: 	Maximum number of iterations. Iteration stops when the iteration count reaches this value.
skip_edge:              Skip edge regions of the survey area when generating CTP.
threshold_Termination:	Iteration convergence threshold. Iteration stops when the difference between current and previous indicators is less than this value.
num_inertia:			Consecutive convergence count. Iteration stops only when convergence signals are triggered consecutively for the given number of times.
'''
CFTM_args = {
    "gird_size_1": 32,
    "gird_size_2": 64,
    "scoring": [0, 0, 0, 0, 0, 0, 2, 4, 8, 16],
    "num_nominated": 10,
    "num_selected": 5,
    "num_max_iterations": 20,
    "pool_size": 16,
    "criterias_1": [2, 3],
    "criterias_2": [1, 3],
    "ratio_inlier_init": 0.25,
    "confidence": 0.9,
    "threshold_RepError": 0.5,
    "threshold_Distance": 3,
    "skip_edge": False,
    "threshold_Termination": 0.001,
    "num_inertia": 3,
}

# Make sure the coordinates of the chunk are in the projected coordinate system;
# here is just an example, this should change to the EPSG code of your projected coordinate system.
reference_args = {
    "PJCS_EPSG": 32647,
    "camera_accuracy_m": 3.0,
    "camera_accuracy_deg": 10.0,
    "marker_accuracy_m": 20.0,
    "marker_accuracy_pix": 0.5,
    "tiepoint_accuracy_pix": 1,
    "camera_reference_enabled": True,
    "marker_reference_enabled": False,
    "rotation_system": 'opk'
}

# 如果想为每期像片的相机参数
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

# Photogrammetry arguments
photogrammetry_args = {
    # Align photos
    "quality_SPC": 1,  # 1(UltraHigh), 2(High), 4(Medium)
    "keypoint_limit": 60000,
    "keypoint_limit_Per_Megapixel": 1000,
    "tiepoint_limit": 0,

    # Error reduction
    "reprojection_error": 0.3,
    "reconstruction_uncertainty": 15,
    "projection_accuracy": 5,
    "max_ratio_remove_points": 0.3,

    # Dense matching
    "quality_DPC": 2,  # 1(UltraHigh), 2(High), 4(Medium)
    "dense_cloud_FileterMode": 'Mild',

    # Build tiled model
    "tiled_model_TiledModelFormat": 'Cesium',
    "tiled_model_pixel_size": 0.3,  # metres
    "tiled_model_face_count": 20000,
    "tiled_model_ModelFormat": 'None',
    "tiled_model_ImageFormat": 'JPEG',
}
# Packing arguments
args = {
    # The Common Feature Track Matching settings
    "CFTM_args": CFTM_args,

    # The process settings
    "process_args": process_args,

    # The photogrammetry settings
    "photogrammetry_args": photogrammetry_args,
    "reference_args": reference_args,
    "bundle_adjustment_args": bundle_adjustment_args,

    # The default settings
    "epoch_mode": epoch_mode,
    "feature_extraction_mode": feature_extraction_mode,
    "begin_from_breakpoint": begin_from_breakpoint,
    "for_all_grids": for_all_grids,
    "reuse_ICTPs": reuse_ICTPs,
    "coalign_chunk_name": coalign_chunk_name,

    # The default path (Recommended to keep the default)
    "code_path": PATH_CODE,
    "project_path": os.path.join(workspace_path, project_name),
    "workspace_path": workspace_path,
    "project_name": project_name,
    "images_path": images_path,
    "check_points_path": check_points_path,
    "unstable_areas_mask_path": unstable_areas_mask_path,
    "process_path": os.path.join(workspace_path, "cftm"),
    "data_package_path": os.path.join(os.path.join(workspace_path, "cftm"), "ProcessData"),
    "feature_path": os.path.join(os.path.join(workspace_path, "cftm"), "FeatureData"),
    "CTP_path": os.path.join(os.path.join(workspace_path, "cftm"), "CommonTiePoints"),
    "iterations_path": os.path.join(os.path.join(workspace_path, "cftm"), "Iterations"),
    "report_path": os.path.join(os.path.join(workspace_path, "cftm"), "Reports"),
    "state_path": os.path.join(os.path.join(workspace_path, "cftm"), "state.json"),
    "arguments_path": os.path.join(os.path.join(workspace_path, "cftm"), "Arguments.json"),
    "state_CFTM_path": os.path.join(os.path.join(workspace_path, "cftm"), "StateCFTM")}
###################################################   END SETUP   ######################################################
########################################################################################################################

# Export and store data from metashape
if __name__ == '__main__':
    starttime0 = time.perf_counter()
    # [0]  Initializing
    # packing arguments
    print(args["process_path"])
    os.makedirs(args["process_path"], exist_ok=True)
    Func_Files.writeArguments(args["arguments_path"], args)
    # load Metashape document and chunk
    state = Func_Breakpoint.check_state(args)
    doc = FuncMs_Basic.openDocument(args)
    chunk = FuncMs_Basic.accessMaxChunk(doc, coalign_chunk_name)
    # Checking compatibility
    found_major_version = ".".join(Metashape.app.version.split('.')[:2])
    if found_major_version != compatible_major_version:
        raise Exception(
            "Incompatible Metashape version: {} != {}".format(found_major_version, compatible_major_version))
    print(state)

    # [1]  Co-alignment
    save_project_each_step = process_args["save_project_each_step"]
    if process_args["process_coalign"] and state['coalign'] == 'NO':
        FuncMs_Basic.alignImages(chunk, quality=photogrammetry_args["quality_SPC"],
                                 KeypointLimit=photogrammetry_args["keypoint_limit"],
                                 KeypointLimitPerMegapix=photogrammetry_args["keypoint_limit_Per_Megapixel"],
                                 TiepointLimit=photogrammetry_args["tiepoint_limit"], realign=True)
        if save_project_each_step:
            doc.save()
    state['coalign'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [2]  CFTM part 1: creat task
    if process_args["process_CTPs_generation"] and state['CFTM_creat_task'] == 'NO':
        Process_CTPsGeneration_part1.creatTask(chunk, args)
        if save_project_each_step:
            doc.save()
    state['CFTM_creat_task'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [3]  Feature Extraction
    if process_args["process_feature_extraction"] and state['feature_extraction'] == 'NO':
        Func_Features.extractFeature(args)
    state['feature_extraction'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [4]  CFTM part 2: run task
    if process_args["process_CTPs_generation"] and state['CFTM_run_task'] == 'NO':
        arguments_path = args["arguments_path"]
        os.system(f"metashape.exe -r {PATH_CODE}/toolbox/Process_CFTM_part2.py {arguments_path}")
    state['CFTM_run_task'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [5]  Iterative Optimization part 1: merge CTPs
    if process_args["process_iterative_optimization"] and state['IO_mergeCTP'] == 'NO':
        Process_IterativeOptimization_part1.mergeCTPs(chunk, args)
        doc.save()
    state['IO_mergeCTP'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [6]  Iterative Optimization part 2: run iterations
    if process_args["process_iterative_optimization"] and state['IO_iterater'] == 'NO':
        Process_IterativeOptimization_part2.iterater(chunk, args)
    state['IO_iterater'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [7]  Iterative Optimization part 3: implement CTPs
    if process_args["process_iterative_optimization"] and state['IO_implementCTPs'] == 'NO':
        chunk = Process_IterativeOptimization_part3.optimizer(chunk, args)
        doc.save()
    state['IO_implementCTPs'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [8]  Analysis CTPs
    if process_args["process_analysis"] and state['analysis_CTPs'] == 'NO':
        Process_Analysis.analysisCTPs(chunk, args)
    state['analysis_CTPs'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [9]  Analysis CPs
    if process_args["process_analysis"] and check_points_path and state['analysis_CPs'] == 'NO':
        Process_Analysis.analysisCPs(chunk, args)
    state['analysis_CPs'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [10]  Analysis Matches
    if process_args["process_analysis_matches"] and state['analysis_matches'] == 'NO':
        Process_Analysis.analysisMatches(args)
    state['analysis_matches'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [11]  Compare Matches
    if process_args["process_compare_matches"] and state['compare_matches'] == 'NO':
        Process_Analysis.compareMatches(chunk, args)
    state['compare_matches'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    # [12]  Reconstruction
    if process_args["process_reconstruction"] and state['reconstruction'] == 'NO':
        Process_Coalign.reconstruction(doc, args)
    state['reconstruction'] = 'DONE'
    Func_Files.writeArguments(args["state_path"], state)

    print('[Script][TimeCost]    All process:', time.perf_counter() - starttime0)
