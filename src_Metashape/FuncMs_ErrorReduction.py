import Metashape
import src_Metashape.FuncMs_Basic as FcBasic


def ErrorReduction(chunk, CriteriasThreshold, bundle_adjustment_args, update=True):
    '''
    input:
        Metashape.document.chunk
        CriteriasThreshold =[0.3, 15, 5], refers to RE,RU,PA respectively
    output:
        Metashape.document.chunk
        removedPointIdList=[point_id,...]
    note:
        The removal will not take effect immediately until optimization is performed or the project file is saved
    '''
    print('[Script]    number of points before Error Reduction:', len(chunk.point_cloud.points))
    points = chunk.point_cloud.points
    CriteriasThreshold_RE = CriteriasThreshold[0]
    CriteriasThreshold_RU = CriteriasThreshold[1]
    CriteriasThreshold_PA = CriteriasThreshold[2]
    MaxRemovePoints = CriteriasThreshold[3]

    # Initialize point quality
    # RE, RU, PA, IC
    MF = Metashape.PointCloud.Filter()
    MF.init(chunk, Metashape.PointCloud.Filter.ReprojectionError)
    PointsCriterion_RE = MF.values
    MF.init(chunk, Metashape.PointCloud.Filter.ReconstructionUncertainty)
    PointsCriterion_RU = MF.values
    MF.init(chunk, Metashape.PointCloud.Filter.ProjectionAccuracy)
    PointsCriterion_PA = MF.values
    # MF.init(chunk, Metashape.PointCloud.Filter.ImageCount)
    # PointsCriterion_IC = MF.values

    # PointsQuality
    RemovedPointIdList = []
    for point_id, point in enumerate(points):
        if PointsCriterion_RE[point_id] > CriteriasThreshold_RE:
            RemovedPointIdList.append(point_id)
        if PointsCriterion_RU[point_id] > CriteriasThreshold_RU:
            RemovedPointIdList.append(point_id)
        if PointsCriterion_PA[point_id] > CriteriasThreshold_PA:
            RemovedPointIdList.append(point_id)
    RemovedPointIdSet = set(RemovedPointIdList)
    if len(RemovedPointIdSet) <= len(points) * 0.3:
        # remove outlier points
        for point_id in RemovedPointIdSet:
            point = chunk.point_cloud.points[point_id]
            point.valid = False
    else:
        # in case there are too much points selected to be removed
        while len(RemovedPointIdSet) > len(points) * MaxRemovePoints:
            CriteriasThreshold_RE += 0.1
            CriteriasThreshold_RU += 2
            CriteriasThreshold_PA += 1
            RemovedPointIdList = []
            for point_id, point in enumerate(points):
                if PointsCriterion_RE[point_id] > CriteriasThreshold_RE:
                    RemovedPointIdList.append(point_id)
                if PointsCriterion_RU[point_id] > CriteriasThreshold_RU:
                    RemovedPointIdList.append(point_id)
                if PointsCriterion_PA[point_id] > CriteriasThreshold_PA:
                    RemovedPointIdList.append(point_id)
            RemovedPointIdSet = set(RemovedPointIdList)
        # remove outlier points
        for point_id in RemovedPointIdSet:
            point = chunk.point_cloud.points[point_id]
            point.valid = False
        print('[Script]    too much points are selected, upper criterias are set automatically:',
              CriteriasThreshold_RE, CriteriasThreshold_RU, CriteriasThreshold_PA)
    # Bundle Adjustment,Estimate tiepoint covariance
    if update:
        FcBasic.bundleAdjustment(chunk, bundle_adjustment_args)
        print('[Script]    number of points after ErrorReduction and optimizeCameras:', len(chunk.point_cloud.points))
    else:
        print('[Script]    number of points after ErrorReduction before optimizeCameras:', len(chunk.point_cloud.points))
    print('[Script]    number of points after ErrorReduction and optimizeCameras:', len(chunk.point_cloud.points))
    return chunk, list(RemovedPointIdSet)
