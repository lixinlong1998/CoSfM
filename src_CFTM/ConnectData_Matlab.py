import numpy as np
import os
import time
import csv
import scipy
import src_CFTM.Func_Files as Func_Files


def buildQueryFiles(QueriesList, Cameras, path):
    '''
    input:
        QueriesList = [Queries1,Queries2...]
            ,where Queries = {keypoint_id:Query}
            ,where Query = [coordinates,QueryProjections]
            ,where coordinates = np.array(Point_Coord_PJCS) with 3 dim
            ,where QueryProjections = [QueryProj,QueryProj,...,QueryProj]
            ,where QueryProj = [u, v, camera_id]
    output:
        QueryCollection = {camera_id:[QueryProj,...]}
            ,where QueryProj = [u, v, keypoint_id]
        questionnaire = {camera_id:[image_path, QueryPoints_path, Answer_path]}
        ./questionnaire file.txt
        ./QueryPoints_{identifier}.txt
        ./QueryPoints_{identifier}.txt
                    ...
        ./QueryPoints_{identifier}.txt
            ,where .txt have two columns for u and v respectively.
    '''
    starttime = time.perf_counter()
    QueryCollection = {}
    for Queries in QueriesList:
        for keypoint_id, Query in Queries.items():
            ValidProjections = Query[1]
            for projection in ValidProjections:
                u = projection[0]
                v = projection[1]
                camera_id = projection[2]
                if camera_id not in QueryCollection:
                    QueryCollection[camera_id] = [[u, v, keypoint_id]]
                else:
                    QueryCollection[camera_id].append([u, v, keypoint_id])

    questionnaire = {}
    # write todolsit
    File = open(path + '/Questionnaire.txt', "w")
    fwriter = csv.writer(File, delimiter=' ', lineterminator='\n')
    QueryPointsFolder = path + '/QueryPoints'
    if not os.path.exists(QueryPointsFolder):
        os.mkdir(QueryPointsFolder)
    for camera_id, QueryProjections in QueryCollection.items():
        image_path = Cameras[camera_id][-2]
        identifier = Cameras[camera_id][-1]
        # 将询问点的文件输出到该文件夹下
        QueryPoints_path = QueryPointsFolder + '/QueryPoints_{0}.txt'.format(identifier)
        Answer_path = QueryPointsFolder + '/Answer_{0}.txt'.format(identifier)

        Func_Files.writeList2D(QueryPoints_path, QueryProjections, cols=[0, 1], Delimiter=' ')
        fwriter.writerow([image_path, QueryPoints_path, Answer_path])
        questionnaire[camera_id] = [image_path, QueryPoints_path, Answer_path]
    File.close()
    print('[Script]        build Query Files[Script][TimeCost]    :', time.perf_counter() - starttime)
    return QueryCollection, questionnaire


def connectImageDescriptor(DPCFeatures1, DPCFeatures2, QueryCollection, questionnaire, camera_ids, CamerasEpoch):
    for image_path, QueryPoints_path, Answer_path in questionnaire:
        identifier = int(Answer_path.split('_')[-1][:-4])
        camera_id = camera_ids[identifier]
        epoch_id = CamerasEpoch[camera_id]
        QueryProjections = QueryCollection[camera_id]
        Answer = Func_Files.readAnswer(Answer_path)
        if epoch_id == 0:
            for Query_id, ImgDes in Answer.items():
                keypoint_id = QueryProjections[Query_id][-1]
                DPCFeatures1[keypoint_id][-1].append(ImgDes)
        elif epoch_id == 1:
            for Query_id, ImgDes in Answer.items():
                keypoint_id = QueryProjections[Query_id][-1]
                DPCFeatures2[keypoint_id][-1].append(ImgDes)
    return DPCFeatures1, DPCFeatures2
