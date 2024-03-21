import csv
import json
import numpy as np

'''
write result to file
'''


def getVariableName(variable):
    loc = locals()
    for k, v in loc.items():
        if loc[k] is variable:
            return k


def writeFile(path):
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter=' ', lineterminator='\n')
    return fwriter


def writeListStatistic(fwriter, VariableName, Variable):
    fwriter.writerow([VariableName + '_AVG', 'MIN', '25P', '50P', '75P', 'MAX'])
    fwriter.writerow([
        '{0:0.5f}'.format(Variable[0]),
        '{0:0.5f}'.format(Variable[1]),
        '{0:0.5f}'.format(Variable[2]),
        '{0:0.5f}'.format(Variable[3]),
        '{0:0.5f}'.format(Variable[4]),
        '{0:0.5f}'.format(Variable[5])])


def writeList2D(path, List2D, cols=[], Delimiter='\t'):
    '''
    input:
        List2D = [list,list,...,list]
            ,where list = [a,b,...,c]
    output:
        .../List.txt
    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter=Delimiter, lineterminator='\n')
    if cols:
        for List in List2D:
            fwriter.writerow([List[i] for i in cols])
    else:
        for List in List2D:
            fwriter.writerow(List)
    File.close()


def readList2D(path, Delimiter='\t'):
    csv.field_size_limit(100000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter=Delimiter, lineterminator='\n')
        List2D = []
        for row in freader:
            List2D.append([eval(i) for i in row])
    return List2D


def writeStatefile(path):
    with open(path, 'w') as file:
        file.write("DONE")


########################################################################################################################
def writeDictionary_Int_List(path, DictItems):
    '''
    input:
        DictItems = {key:value} , where key is int and value is list
        value = [list, np.array, np.ndarray]
    output:
        .../DictItems.txt

    适用于此函数的变量：
        Cameras = {camera_id:Camera}
            ,where Camera = [transform, reference, covariance, sensors_id, state, path, identifier]
            ,where transform = [camera_transform, camera_loc, camera_rot, camera_transform_PJCS]
            ,where reference = [camera_loc_Ref,location_accuracy,camera_rot_Ref,rotation_accuracy]
            ,where covariance = [camera_location_covariance, camera_rotation_covariance]
            ,where states = [label,camera_loc_Ref,camera_rot_Ref, selected, camera_orientation]
        Sensors = {sensor_id:Sensor}
            ,where calibration_params = [f, cx, cy, b1, b2,  k1, k2, k3, k4,  p1, p2, p3, p4,  height, width]
            ,where calibration_matrix = nxn covariance matrix
        DPCFeatures = {keypoint_id:DPCFeature}
            ,where feature = [coordinates,geoDescriptor,imgDescriptor]
            ,where coordinates = np.array(Point_Coord_LCCS) with 3 dim
            ,where geoDescriptor = np.array(geoDescriptor) with 33 dim
            ,where imgDescriptor = []
        Queries = {keypoint_id:Query}
            ,where Query = [coordinates,QueryProjections]
            ,where coordinates = np.array(Point_Coord_PJCS) with 3 dim
            ,where QueryProjections = [] it is empty here, but should be [QueryProj,QueryProj,...,QueryProj]
            ,where QueryProj = [u, v, camera_id]
        QueryCollection = {camera_id:[QueryProj,...]}
            ,where QueryProj = [u, v, keypoint_id]
    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    for int_id, value_list in DictItems.items():
        fwriter.writerow([int_id, value_list])
    File.close()


def readDictionary_Int_List(path):
    csv.field_size_limit(100000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        DictItems = {}
        for row in freader:
            int_id = int(row[0])
            value_list = eval(row[1].replace('array', 'np.array').replace('Matrix', 'np.matrix'))
            DictItems[int_id] = value_list
    return DictItems


def writeDictionary_Int_Int(path, DictItems):
    '''
    input:
        DictItems = {key:value} ,the key is int while the value is also the int
        value = [list, np.array, np.ndarray]
    output:
        .../DictItems.txt

    适用于此函数的变量：
        camera_ids = {identifier:camera_id}

    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    for int_id, int_value in DictItems.items():
        fwriter.writerow([int_id, int_value])
    File.close()


def readDictionary_Int_Int(path):
    csv.field_size_limit(100000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        DictItems = {}
        for row in freader:
            int_id = int(row[0])
            int_value = int(row[1])
            DictItems[int_id] = int_value
    return DictItems


def writeDictionary_Tuple_List(path, DictItems):
    '''
    input:
        DictItems = {key:value} ,the key is tuple while the value is list
        value = [list, np.array, np.ndarray]
    output:
        .../DictItems.txt

    适用于此函数的变量：
        Cubes = {grid_id:conners}
            ,where conners = [topLeft, topRight, leftBottom, rightBottom, upper, lower]
        Camera_IdList_cube = {grid_id:Camera_IdList}
            ,where Camera_IdList = [camera_id,camera_id,...,camera_id]
        MarkersGrid = {(r, c):[Marker,...]}
            ,where Marker = [MarkerPoint,track,Information]
    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    for tuple_id, value_list in DictItems.items():
        fwriter.writerow([tuple_id, value_list])
    File.close()


def readDictionary_Tuple_List(path):
    csv.field_size_limit(100000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        DictItems = {}
        for row in freader:
            tuple_id = eval(row[0])
            value_list = eval(row[1].replace('array', 'np.array').replace('Matrix', 'np.matrix'))
            DictItems[tuple_id] = value_list
    return DictItems


def writeDictionary_Str_Int(path, DictItems):
    '''
    input:
        DictItems = {key:value} ,the key is string while the value is int
        value = [list, np.array, np.ndarray]
    output:
        .../DictItems.txt

    适用于此函数的变量：
        cameraPaths = {camera_photo_path:identifier}

    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    for key_str, value_int in DictItems.items():
        fwriter.writerow([key_str, value_int])
    File.close()


def readDictionary_Str_Int(path):
    csv.field_size_limit(100000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        DictItems = {}
        for row in freader:
            key_str = str(row[0])
            value_int = int(row[1])
            DictItems[key_str] = value_int
    return DictItems


########################################################################################################################
def writeArguments(path, Dictionary):
    '''
    input:
        Dictionary = {'VariableName':variable,...}
    output:
        .../Arguments.txt
    '''
    with open(path, 'w') as json_file:
        json.dump(Dictionary, json_file)


def readArguments(path):
    with open(path, 'r') as json_file:
        Dictionary = json.load(json_file)
    return Dictionary


def writeCommonTracks(path, CommonTracks, CommonTracksMatches):
    '''
    input:
        CommonTracks = [CommonTrack,CommonTrack,...,CommonTrack]
            ,where CommonTrack = [view,view,...,view]
            ,where view = [camera_id, projection_info, index_info]
            ,where projection_info = [u, v]
            ,where index_info = [identifier, keypoints_id, originTrack_id]
        CommonTracksMatches = [CommonTrackMatches,CommonTrackMatches,...,CommonTrackMatches]
            ,where CommonTrackMatches = [CrossMatch,…]
            ,where CrossMatch = [feature1, feature2, MatchQuality]
            ,where feature = [camera_id, keypoint_id]
            ,where MatchQuality = [distance,angles],distance is the L2 normal distance between descriptors of match.
    output:
        .../ICTPsTracks.txt
    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    for i, CommonTrack in enumerate(CommonTracks):
        fwriter.writerow([CommonTrack, CommonTracksMatches[i]])
    File.close()


def readCommonTracks(path):
    csv.field_size_limit(100000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        CommonTracks = []
        CommonTracksMatches = []
        for row in freader:
            CommonTracks.append(eval(row[0]))
            CommonTracksMatches.append(eval(row[1]))
    return CommonTracks, CommonTracksMatches


def writeCamerasEpoch(path, CamerasEpoch):
    '''
    input:
        CamerasEpoch = [epoch_id,epoch_id,...,epoch_id] indexed by camera_id
    output:
        .../CamerasEpoch.txt
    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    for camera_id, epoch_id in enumerate(CamerasEpoch):
        fwriter.writerow([camera_id, epoch_id])
    File.close()


def readCamerasEpoch(path):
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        CamerasEpoch = []
        for row in freader:
            CamerasEpoch.append(int(row[1]))
    return CamerasEpoch


def writeCoordinate(path, CoordinateTransform, CoordinateAttribute):
    '''
    input:
        CoordinateTransform = [M, LSE, T, R, S,M_inv]
        CoordinateAttribute = [crs_name, crs_PJCS_wkt, crs_GCCS_wkt, crs_GDCS_wkt]
    output:
        .../Coordinate.txt
    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    fwriter.writerow([CoordinateTransform])
    fwriter.writerow([CoordinateAttribute])
    File.close()


def readCoordinate(path):
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        Coordinate = []
        for row in freader:
            Coordinate.append(row[0].replace('array', 'np.array').replace('Matrix', 'np.matrix'))
        CoordinateTransform = eval(Coordinate[0])
        CoordinateAttribute = eval(Coordinate[1])
    return CoordinateTransform, CoordinateAttribute


def writeMatchesReport(path, MatchesReport_cube):
    '''
    input:
        MatchesReport_cube = [pair,pair,...,pair]
            ,where pair = [(camera_id1, camera_id2), Matches, inlierMatches]
            ,where Matches = np.array([[keypoint_id, keypoint_id],…]).shape(matchesNum,2)
            ,where inlierMatches = np.array([[keypoint_id, keypoint_id],…]).shape(matchesNum,2)
    output:
        .../MatchesReport_.txt
    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    for pair in MatchesReport_cube:
        fwriter.writerow([pair])
    File.close()


def readMatchesReport(path):
    '''
    input:
        .../MatchesReport_.txt
    output:
        MatchesReport_cube = [pair,pair,...,pair]
            ,where pair = [(camera_id1, camera_id2), Matches, inlierMatches]
            ,where Matches = np.array([[keypoint_id, keypoint_id],…]).shape(matchesNum,2)
            ,where inlierMatches = np.array([[keypoint_id, keypoint_id],…]).shape(matchesNum,2)
    '''
    csv.field_size_limit(100000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        MatchesReport_cube = []
        for row in freader:
            MatchesReport_cube.append(eval(row[0].replace('array', 'np.array').replace('Matrix', 'np.matrix')))
    return MatchesReport_cube


def writeMarkersMaskbyUnstable(path, MarkersGrid):
    '''
    input:
        MarkersGrid = {(r, c):[Marker,...]}
            ,where Marker = [MarkerPoint,track,Information]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
            ,where track = [view,view,...,view]
            ,where Information = [[point_id,Track_id],MarkerKey] or [[CrossMatch,CrossMatch,...,CrossMatch],MarkerKey]
    output:
        .../MarkersGrid_MaskbyUnstable.txt
    '''
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    # fwriter.writerow(['X1', 'Y1'])
    for grid_id, Markers in MarkersGrid.items():
        for i, Marker in enumerate(Markers):
            fwriter.writerow(['{0:0.5f}'.format(Marker[0][1][0]), '{0:0.5f}'.format(Marker[0][1][1])])
    File.close()


def readMarkersMaskbyUnstable(path):
    '''
    input:
        .../MarkersGrid_MaskbyUnstable.txt
    output:
        MarkersGrid = {(r, c):[Marker,...]}
            ,where Marker = [MarkerPoint,track,Information]
            ,where MarkerPoint = [Point_Coord_CLCS,Point_Coord_PJCS]
            ,where track = [view,view,...,view]
            ,where Information = [[point_id,Track_id],MarkerKey] or [[CrossMatch,CrossMatch,...,CrossMatch],MarkerKey]
    '''
    csv.field_size_limit(100000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        MarkersGrid_MaskbyUnstable = []
        for row in freader:
            MarkersGrid_MaskbyUnstable.append(eval(row[0].replace('array', 'np.array').replace('Matrix', 'np.matrix')))
    return MarkersGrid_MaskbyUnstable


def creatIterationsMERE(path):
    File = open(path, "w")
    fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
    fwriter.writerow(["Iteration state: running"])
    File.close()


def readIterationsMERE(path):
    with open(path, "r") as file:
        freader = csv.reader(file, delimiter='\t', lineterminator='\n')
        headers = next(freader)
        record_MERE = {}
        for row in freader:
            record_MERE[row[0]] = row[1]
    return record_MERE


########################################################################################################################


def readQuestionnaire(path):
    csv.field_size_limit(100000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter=' ', lineterminator='\n')
        questionnaire = []
        for row in freader:
            image_path = str(row[0])
            QueryPoints_path = str(row[1])
            Answer_path = str(row[2])
            questionnaire.append([image_path, QueryPoints_path, Answer_path])
    return questionnaire


def readAnswer(path):
    csv.field_size_limit(100000000)
    with open(path) as f:
        freader = csv.reader(f, delimiter='\t', lineterminator='\n')
        Answer = {}
        for row in freader:
            Query_id = int(row[0])
            Answer[Query_id] = eval(row[1])
    return Answer

# def writeKeypoints(path, chunkKeypoints):
#     '''
#     input:
#         chunkKeypoints = {identifier:2darray.shape(numKP,6).dtype(float32),...}
#
#     output:
#         .../Keypoints.txt
#     '''
#     File = open(path, "w")
#     fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
#     for identifier, Keypoints_camera in chunkKeypoints.items():
#         fwriter.writerow([identifier, list(Keypoints_camera)])
#     File.close()
#
#
# def readKeypoints(path):
#     with open(path) as f:
#         freader = csv.reader(f, delimiter='\t', lineterminator='\n')
#         CoordinateTransform = eval(freader[0])
#         CoordinateAttribute = eval(freader[1])
#     return chunkKeypoints
#
#
# def writeDescriptors(path, chunkDescriptors):
#     '''
#     input:
#         chunkDescriptors = {identifier:2darray.shape(numKP,128).dtype(uint8),...}
#     output:
#         .../Descriptors.txt
#     '''
#     File = open(path, "w")
#     fwriter = csv.writer(File, delimiter='\t', lineterminator='\n')
#     for identifier, Descriptors_camera in chunkDescriptors.items():
#         fwriter.writerow([identifier, list(Descriptors_camera)])
#     File.close()
#
#
# def readDescriptors(path):
#     with open(path) as f:
#         freader = csv.reader(f, delimiter='\t', lineterminator='\n')
#         CoordinateTransform = eval(freader[0])
#         CoordinateAttribute = eval(freader[1])
#     return chunkDescriptors

#
# def writeDatabase(path, *args):
#     conn = sqlite3.connect(path)  # datapackage.db
#     cursor = conn.cursor()
#     cursor.execute('''CREATE TABLE Cubes(id INTEGER PRIMARY KEY, grid_id TEXT, data BLOB)''')
#     cursor.execute('''CREATE TABLE Cube_CameraIdList
#                       (id INTEGER PRIMARY KEY, name TEXT, data BLOB)''')
#     cursor.execute('''CREATE TABLE Cameras
#                       (id INTEGER PRIMARY KEY, name TEXT, data BLOB)''')
#     cursor.execute('''CREATE TABLE camera_ids
#                       (id INTEGER PRIMARY KEY, name TEXT, data BLOB)''')
#     cursor.execute('''CREATE TABLE Sensors
#                       (id INTEGER PRIMARY KEY, name TEXT, data BLOB)''')
#     cursor.execute('''CREATE TABLE Coordinate
#                       (id INTEGER PRIMARY KEY, name TEXT, data BLOB)''')
#     cursor.execute('''CREATE TABLE Keypoints
#                       (id INTEGER PRIMARY KEY, name TEXT, data BLOB)''')
#     cursor.execute('''CREATE TABLE Descriptors
#                       (id INTEGER PRIMARY KEY, name TEXT, data BLOB)''')
#     conn.commit()
