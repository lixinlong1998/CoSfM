import os
def splitChunk_individual(coalign_chunk, epochs, camera_ids, CamerasEpoch, epoch_mode):
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
                if os.path.basename(os.path.dirname(camera_path)) != epochs[i]:
                    camerasRemoved_list.append(camera)

        chunk.remove(camerasRemoved_list)


def splitChunk_pairwise(coalign_chunk, chunksName, epochs, camera_ids, CamerasEpoch, epoch_mode):
    '''
    function:
        split the all co-aligned chunk into pairwise co-aligned chunks.
        note: this function only used for case of all co-aligned with triple epochs
    parameters should be given as follows:
        coalign_chunk
        epochs = [['2020:10:14'],['2021:07:03'],['2021:12:03']]
        chunksName = ['coalign_12', 'coalign_13', 'coalign_23']
    return:
        coalign_12, coalign_13, coalign_23
    '''
    # duplicate chunk
    chunksNum = len(epochs)
    chunksList = []
    for i in range(chunksNum):
        chunkName = chunksName[i]
        locals()[chunkName] = coalign_chunk.copy(keypoints=False)
        locals()[chunkName].label = chunksName[i]
        chunksList.append(locals()[chunkName])
    # remove cameras which capture date is not in line with the epoch of chunk
    for i, chunk in enumerate(chunksList):
        # the epoch that need removed
        epoch = epochs[-1 * (i + 1)]
        # remove cameras from epoch above
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
                if os.path.basename(os.path.dirname(camera_path)) != epochs[i]:
                    camerasRemoved_list.append(camera)
        chunk.remove(camerasRemoved_list)
