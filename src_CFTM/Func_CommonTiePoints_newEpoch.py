def getEpoch_byDATE(camera_ids):
    '''
    This function can automatically classify the epoch of images.
    This function will search the date of images firstly and find out how many epochs are there,
    it returns:
        epochs = [['2020:10:14'], ['2021:07:03'], ['2021:12:03'], ['2021:12:06']]
    Note that this function define images took within one day as an epoch image data sets, those image data set acquired
    in several days are not suitable for this function.
    '''
    datelist = []
    for identifier in camera_ids:
        datelist.append(int(str(identifier)[:8]))
    epochlist = list(set(datelist))
    epochlist.sort()
    epochs = []
    for epoch in epochlist:
        epochs.append([str(epoch)[0:4] + ':' + str(epoch)[4:6] + ':' + str(epoch)[6:8]])
    return epochs


def analyseSignal_Epoch(epochNum):
    '''
    Using 0 or 1 to define whether this epoch has images.
        input: epochNum = [e1,e2,e3,e4]
        output: Signal = [0,0,1,1]
    '''
    signal = [0] * len(epochNum)
    for i, ei in enumerate(epochNum):
        if ei > 0:
            signal[i] = 1
    return signal


def analyseCameras_Epoch(camera_ids, epochs):
    '''
    input:
        camera_ids = {identifier:camera_id}
        epochs = [['2020:10:14'], ['2021:07:03'], ['2021:12:03'], ['2021:12:06']]
    output:
        CamerasEpoch = [epoch_id,epoch_id,...,epoch_id] indexed by camera_id
    assign which epochs dose camera belongs to.
    '''
    CamerasEpoch = [-1] * len(camera_ids)
    for identifier, camera_id in camera_ids.items():
        # get the original datetime of image
        identifier_Date = str(identifier)[:8]  # 20210404
        Date = identifier_Date[0:4] + ':' + identifier_Date[4:6] + ':' + identifier_Date[6:8]  # '2021:04:04'
        # check which epoch dose the photo belongs to
        for i, epoch in enumerate(epochs):
            if Date in epoch:
                CamerasEpoch[camera_id] = i
    if -1 in CamerasEpoch:
        raise Exception("[Script]    the given epochs are not include all dates")
    return CamerasEpoch


def analyseTracks_Epoch(Tracks, epochs, CamerasEpoch):
    '''
    input:
        Tracks = [Track, Track,…, Track] (index by track_id)
            ,where Track = [view,view,...view],
            ,where view = [camera_id,projection_info],
            ,where projection_info = [projection.coord[0], projection.coord[1], projection.size, projection.track_id].
        epochs = [['2020:10:14'], ['2021:07:03'], ['2021:12:03'], ['2021:12:06']]
            each item is a list contain the date of photographing,it could contains more than one date in an epoch.
    output:
        TracksEpoch = [[e1,e2,e3,e4], [e1,e2,e3,e4],…, [e1,e2,e3,e4]]
            ,where e1 means the number of views taken from epoch1 in a Track ,and as well as e2,e3,e4...
    This function could analysis the number of images for each epoch in a Track.
    '''
    TracksEpoch = []
    for Track in Tracks:
        epochNum = [0] * len(epochs)
        for camera_id, projection_info in Track:
            epochNum[CamerasEpoch[camera_id]] += 1
        TracksEpoch.append(epochNum)
    return TracksEpoch


def getICTPsIndex(Points, point_ids, TracksEpoch):
    '''
    input:
        Points={point_id:Point,...}
            ,where Point = [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, track_id]
            Coordinates=[[X, Y, Z],[X, Y, Z]]   # in CLCS & PJCS respectively
            Covariance = 3x3 covariance matrix
            StandardDeviation = [stdX, stdY, stdZ]  # in mm
            Quality = [RE, RU, PA, IC]
        point_ids = [point_id,point_id,..,point_id] by track_id
        TracksEpoch = [[e1,e2,e3,e4], [e1,e2,e3,e4],…, [e1,e2,e3,e4]] by track_id
    output:
        ICTPs_IdList = [point_id,point_id,...,point_id]
        ICTPs_Signal = [signal,signal,...,signal]
            ,where signal = [0,1,0,1]*len(epochs)
    This function analysis Common Tie Points(CTPs) based on TracksEpoch.
    CTP is defined as the point, which corresponding Track contains images from different epochs
    '''
    ICTPs_IdList = []
    ICTPs_Signal = []
    for track_id, epochNum in enumerate(TracksEpoch):
        point_id = point_ids[track_id]
        # skip invalid Track
        if point_id == -1 or Points[point_id][5] == 0:
            continue
        signal = analyseSignal_Epoch(epochNum)
        # if the Track contains images from more than one epoch, this Track is considered as EpochTrack
        if sum(signal) > 1:
            ICTPs_IdList.append(point_id)
            ICTPs_Signal.append(signal)
    return ICTPs_IdList, ICTPs_Signal


def getICTPsTracks(Tracks, Points, epochs, CamerasEpoch, ICTPs_IdList):
    '''
    input:
        Tracks = [Track, Track,…, Track] (index by track_id)
            ,where Track = [view,view,...view],
            ,where view = [camera_id,projection_info],
            ,where projection_info = [projection.coord[0], projection.coord[1], projection.size, projection.track_id]
        Points={point_id:Point,...}
            ,where Point = [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, track_id]
            Coordinates=[[X, Y, Z],[X, Y, Z]]   # in CLCS & PJCS respectively
            Covariance = 3x3 covariance matrix
            StandardDeviation = [stdX, stdY, stdZ]  # in mm
            Quality = [RE, RU, PA, IC]
        epochs = [['2020:10:14'], ['2021:07:03'], ['2021:12:03'], ['2021:12:06']]
        CamerasEpoch = [epoch_id,epoch_id,...,epoch_id] indexed by camera_id
        ICTPs_IdList = [point_id,point_id,...,point_id]
    output:
        ICTPsTracks = [EpochTrack, EpochTrack,…, EpochTrack] (with index accompanied by ICTPsTracks_ids)
            ,where EpochTrack = [views_e1,views_e2,views_e3,views_e4],
            ,where views_e1 = [view,view,...,view]
            ,where view = [camera_id,projection_info],
            ,where projection_info = [projection.coord[0], projection.coord[1], projection.size, projection.track_id].
        ICTPsTracks_ids = [track_id,track_id,...,track_id]
    This function analysis epoch Track(ETs) based on TracksEpoch.
    ETs is defined as the Track that contains images from different epochs
    '''
    ICTPsTracks = []
    ICTPsTracks_ids = []
    for point_id in ICTPs_IdList:
        track_id = Points[point_id][-1]
        Track = Tracks[track_id]
        EpochTrack = [[] for i in range(len(epochs))]
        for view in Track:
            camera_id = view[0]
            EpochTrack[CamerasEpoch[camera_id]].append(view)
        ICTPsTracks.append(EpochTrack)
        ICTPsTracks_ids.append(track_id)
    return ICTPsTracks, ICTPsTracks_ids


def getICTPsCoord(Points, ICTPs_IdList):
    '''
    get the coordinates of Common Tie Points(CTPs)
    parameters:
        Points={point_id:Point,...}
            ,where Point = [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, track_id]
            Coordinates=[[X, Y, Z],[X, Y, Z]]   # in CLCS & PJCS respectively
            Covariance = 3x3 covariance matrix
            StandardDeviation = [stdX, stdY, stdZ]  # in mm
            Quality = [RE, RU, PA, IC]
        ICTPs_IdList = [point_id,...]
    return:
        ICTPsCoord_CLCS = [[X,Y,Z],[X,Y,Z],...,[X,Y,Z]]  # in CLCS
        ICTPsCoord_PJCS = [[X,Y,Z],[X,Y,Z],...,[X,Y,Z]]  # in PJCS
    '''
    ICTPsCoord_CLCS = []
    ICTPsCoord_PJCS = []
    for point_id in ICTPs_IdList:
        Point = Points[point_id]
        Coordinates_CLCS = Point[0][0]
        Coordinates_PJCS = Point[0][1]
        ICTPsCoord_CLCS.append([Coordinates_CLCS[0], Coordinates_CLCS[1], Coordinates_CLCS[2]])
        ICTPsCoord_PJCS.append([Coordinates_PJCS[0], Coordinates_PJCS[1], Coordinates_PJCS[2]])
    return ICTPsCoord_CLCS, ICTPsCoord_PJCS


def getMatches(chunk, point_ids, Tracks):
    '''
    Pair_ids = {}:
        dict that return [all, valid, invalid] by given (camera_id,camera_id);
    Match_ids = {}:
        dict that return [valid matches, invalid matches] by given (camera_id,camera_id);
        where matches = [match,match,...,match],
        where match = (projection_info,projection_info)
        where projection_info = [projection.coord[0], projection.coord[1], projection.size].

    Note that the Pair_ids and Match_ids is completely the same with the matches viewed in GUI.
    To our knowledge, the matches here are only those used to construct Track, and the valid match means the Track
    triangulated successfully.
    So the valid or invalid flag only refers whether this match is USED, the process of geometric verified, misalign of
    Track and RANSAC based Track filtering might not considered here.
    '''
    cameras = chunk.cameras
    points = chunk.point_cloud.points

    # Creat full connection pairs
    Pair_ids = {}
    Match_ids = {}
    for camera_id1 in range(len(cameras) - 1):
        for camera_id2 in range(camera_id1 + 1, len(cameras)):
            Pair_ids[(camera_id1, camera_id2)] = [0, 0, 0]  # [all, valid, invalid]
            Match_ids[(camera_id1, camera_id2)] = [[], []]  # [valid matches, invalid matches]

    for track_id, Track in enumerate(Tracks):
        point_id = point_ids[track_id]
        for i, view1 in enumerate(Track[:-1]):
            for view2 in Track[i + 1:]:
                camera_id1 = view1[0]
                camera_id2 = view2[0]
                if camera_id1 < camera_id2:
                    # pair = (camera_epoch1,camera_epoch2)
                    pair = (camera_id1, camera_id2)
                    # match = ([proj.coordu,proj.coordv,proj.size],[proj.coordu,proj.coordv,proj.size])
                    match = (view1[1], view2[1])
                elif camera_id1 > camera_id2:
                    pair = (camera_id2, camera_id1)
                    match = (view2[1], view1[1])
                else:
                    raise Exception("[Script]    creating pair:a Track contains same camera!")

                # add matches number to dirt
                Pair_ids[pair][0] += 1
                if point_id == -1 or not points[point_id].valid:
                    Pair_ids[pair][2] += 1
                    Match_ids[pair][1].append(match)
                else:
                    Pair_ids[pair][1] += 1
                    Match_ids[pair][0].append(match)
    return Pair_ids, Match_ids


def getCommonMatches(epochPair, cameras, points, point_ids, CamerasEpoch, ICTPsTracks, ICTPsTracks_ids):
    '''
    parameters:
        epochPair = [1,2]   specify which pair needs to assess ,where 0,1,2 indicate the first, second, third epoch
    return:
        EpochPairs = {}:
            dict that return [all, valid, invalid] by given (camera_id,camera_id);
        CommonMatches = {}:
            dict that return [valid matches, invalid matches] by given (camera_id,camera_id);
            where matches = [match,match,...,match],
            where match = (projection_info,projection_info)
            where projection_info = [projection.coord[0], projection.coord[1], projection.size, projection.track_id].

    Note that the valid or invalid flag only refers whether this match is USED, the process of geometric verified,
    misalign of Track and RANSAC based Track filtering might not considered here.
    '''
    epoch_first, epoch_second = epochPair

    # Creat full connection pairs
    EpochPairs = {}
    CommonMatches = {}
    for camera_id1, camera1 in enumerate(cameras):
        if CamerasEpoch[camera_id1] == epoch_first:
            for camera_id2, camera2 in enumerate(cameras):
                if CamerasEpoch[camera_id2] == epoch_second:
                    EpochPairs[(camera_id1, camera_id2)] = [0, 0, 0]  # [all, valid, invalid]
                    CommonMatches[(camera_id1, camera_id2)] = [[], []]  # [valid matches, invalid matches]

    # analysis matches from epoch Tracks
    for i, Track in enumerate(ICTPsTracks):
        track_id = ICTPsTracks_ids[i]
        point_id = point_ids[track_id]
        for i, view1 in enumerate(Track[:-1]):
            for view2 in Track[i + 1:]:
                camera_id1 = view1[0]
                camera_id2 = view2[0]
                if camera_id1 < camera_id2:
                    # pair = (camera_epoch1,camera_epoch2)
                    pair = (camera_id1, camera_id2)
                    # match = ([proj.coordu,proj.coordv,proj.size],[proj.coordu,proj.coordv,proj.size])
                    match = (view1[1], view2[1])
                elif camera_id1 > camera_id2:
                    pair = (camera_id2, camera_id1)
                    match = (view2[1], view1[1])
                else:
                    # skip pairs not across epochs
                    continue
                # add matches number to dirt
                EpochPairs[pair][0] += 1
                if point_id == -1 or not points[point_id].valid:
                    EpochPairs[pair][2] += 1
                    CommonMatches[pair][1].append(match)
                else:
                    EpochPairs[pair][1] += 1
                    CommonMatches[pair][0].append(match)
    return EpochPairs, CommonMatches
