from config import *
import Metashape
from src_CFTM import ConnectData_Metashape
from src_CFTM import Func_CommonTiePoints

'''Please run this script using the following command in cmd：

metashape.exe -r <...\CFTM\src_Process\Process_SplitChunk.py>

Introduction:
确定coalign chunk（指定chunk_name，或者程序默认将像片数最大的chunk作为coalign chunk）
根据像片的拍摄日期将coalign中像片按‘天’分组（也可以手动指定epochs），每组像片根据其日期先后顺序划分到‘chunk_i’
# epochs = [['2020:10:14'], ['2021:07:03'], ['2021:12:03'], ['2021:12:06']]
'''
#################################################       SETUP      #####################################################
project_path = r"K:\Baige\Aprojects\Coalign_e3e5\Baige_e3e5_Coalign_Reconstruction.psx"
chunk_name = ''
epoch_mode = "DATE"  # or "FOLDER"


#################################################   END OF SETUP   #####################################################

def splitChunk(coalign_chunk, epochs, epoch_mode):
    '''
    function:
        split the coalign chunk into individual chunks with its own epoch images
    parameters should be given as follows:
        coalign_chunk
        epochs = [['2020:10:14'],['2021:07:03'],['2021:12:03']]
    return:
        chunk_1,chunk_2,chunk_3
    '''
    # duplicate chunk
    chunksNum = len(epochs)
    chunksList = []
    for i in range(chunksNum):
        chunkName = 'chunk_' + str(i)
        locals()[chunkName] = coalign_chunk.copy(keypoints=False)
        locals()[chunkName].label = 'chunk_' + str(i)
        chunksList.append(locals()[chunkName])
    # remove cameras which capture date is not in line with the epoch of chunk
    for i, chunk in enumerate(chunksList):
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

    # [2]  access chunk
    if chunk_name:
        # choose the chunk with given name
        for chunk_i in Metashape.app.document.chunks:
            if chunk_i.label == chunk_name:
                chunk = chunk_i
    else:
        # choose the chunk with maximum images
        chunk_images_num = 0
        for chunk_i in Metashape.app.document.chunks:
            chunki_images_num = len(chunk_i.cameras)
            if chunki_images_num >= chunk_images_num:
                chunk_images_num = chunki_images_num
                chunk = chunk_i
            else:
                continue

    # [3]  access chunk attributes
    camera_ids = ConnectData_Metashape.getCameraIds(chunk)
    cameraPaths = ConnectData_Metashape.getCameraPaths(chunk)
    if epoch_mode == "DATE":
        epochs = Func_CommonTiePoints.getEpoch_byDATE(camera_ids)
    elif epoch_mode == "FOLDER":
        epochs = Func_CommonTiePoints.getEpoch_byFOLDER(cameraPaths)

    # [4]  split coaligned chunk into individual chunks
    splitChunk(chunk, epochs, epoch_mode)

    # [5] save project
    doc.save()
