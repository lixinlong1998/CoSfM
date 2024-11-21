import shutil
import os
import src_CFTM.Func_Files as Func_Files


def check_state(args):
    state_path = args["state_path"]
    begin_from_breakpoint = args["begin_from_breakpoint"]
    pre_defined_state = {'coalign': 'NO',
                         'CFTM_creat_task': 'NO',
                         'feature_extraction': 'NO',
                         'CFTM_run_task': 'NO',
                         'IO_mergeCTP': 'NO',
                         'IO_iterater': 'NO',
                         'IO_implementCTPs': 'NO',
                         'analysis_CTPs': 'NO',
                         'analysis_CPs': 'NO',
                         'analysis_matches': 'NO',
                         'compare_matches': 'NO',
                         'reconstruction': 'NO',
                         'finished_grids': [],
                         'finished_iterations': []
                         }

    if not os.path.exists(state_path):
        Func_Files.writeArguments(state_path, pre_defined_state)
        print(f"Empty state file '{state_path}' created.")
    state = Func_Files.readArguments(state_path)

    if begin_from_breakpoint:
        return state
    else:
        # clear the content in folder
        if os.path.isdir(args["data_package_path"]):
            shutil.rmtree(args["data_package_path"])
        # if os.path.isdir(args["feature_path"]):
        #     shutil.rmtree(args["feature_path"])
        if os.path.isdir(args["CTP_path"]):
            shutil.rmtree(args["CTP_path"])
        if os.path.isdir(args["iterations_path"]):
            shutil.rmtree(args["iterations_path"])
        if os.path.isdir(args["report_path"]):
            shutil.rmtree(args["report_path"])
        return pre_defined_state
