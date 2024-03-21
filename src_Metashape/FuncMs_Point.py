import Metashape
import math
import src_Metashape.FuncMs_Transform as MsTrans


def getPointsCriterion(chunk, IdList=[]):
    '''
    input:
        Metashape.document.chunk
        optional parameters:IdList is a set of point_id
            Note that the index of returned PointQuality are not align with point_id list but the given IdList.
    output:
        PointsCriterion = [PointCriterion,PointCriterion,...,PointCriterion]
            ,where PointCriterion = [RE, RU, PA, IC]
    optional parameters:
        IdList is a set of point_id
        Note that the index of returned PointQuality are not align with point_id list but the given IdList.
    '''
    points = chunk.point_cloud.points

    # Initialize point quality
    # RE, RU, PA, IC
    MF = Metashape.PointCloud.Filter()
    MF.init(chunk, Metashape.PointCloud.Filter.ReprojectionError)
    PointsCriterion_RE = MF.values
    MF.init(chunk, Metashape.PointCloud.Filter.ReconstructionUncertainty)
    PointsCriterion_RU = MF.values
    MF.init(chunk, Metashape.PointCloud.Filter.ProjectionAccuracy)
    PointsCriterion_PA = MF.values
    MF.init(chunk, Metashape.PointCloud.Filter.ImageCount)
    PointsCriterion_IC = MF.values

    # PointsQuality
    if IdList:
        PointsCriterion = []
        for point_id in IdList:
            PointCriterion = [-1, -1, -1, -1]
            PointCriterion[0] = PointsCriterion_RE[point_id]
            PointCriterion[1] = PointsCriterion_RU[point_id]
            PointCriterion[2] = PointsCriterion_PA[point_id]
            PointCriterion[3] = PointsCriterion_IC[point_id]
            PointsCriterion.append(PointCriterion)
        return PointsCriterion
    else:
        PointsCriterion = []
        for point_id, point in enumerate(points):
            PointCriterion = [-1, -1, -1, -1]
            PointCriterion[0] = PointsCriterion_RE[point_id]
            PointCriterion[1] = PointsCriterion_RU[point_id]
            PointCriterion[2] = PointsCriterion_PA[point_id]
            PointCriterion[3] = PointsCriterion_IC[point_id]
            PointsCriterion.append(PointCriterion)
        return PointsCriterion


def getPointsCov_value(chunk, IdList=[]):
    '''
    input:
        Metashape.document.chunk
        optional parameters:IdList is a set of point_id
            Note that the index of returned PointQuality are not align with point_id list but the given IdList.
    output:
        PointsCov = [PointCov,PointCov,...,PointCov]
            ,where PointCov = [covXX(m2), covXY(m2), covXZ(m2),
                               covYX(m2), covYY(m2), covYZ(m2),
                               covZX(m2), covZY(m2), covZZ(m2)]
    This function must follow the step of 'optimization' with option 'estimating points covariance' checked.
    '''
    points = chunk.point_cloud.points
    M, T, R = MsTrans.getCoordTransMat(chunk)
    # covXX(m2), covXY(m2), covXZ(m2),
    # covYX(m2), covYY(m2), covYZ(m2),
    # covZX(m2), covZY(m2), covZZ(m2).
    PointsCov = [[-1, -1, -1, -1, -1, -1, -1, -1, -1] for i in range(len(points))]
    for point_id, point in enumerate(points):
        # skipping invalid points
        if not point.valid:
            continue
        # Transform point covariance matrix from CLCS to PJCS   # math.sqrt(point_cov[0, 0]) * 1000
        point_cov = R * point.cov * R.t()
        PointsCov[point_id][0] = point_cov[0, 0]  # covXX(m2)
        PointsCov[point_id][1] = point_cov[0, 1]
        PointsCov[point_id][2] = point_cov[0, 2]
        PointsCov[point_id][3] = point_cov[1, 0]
        PointsCov[point_id][4] = point_cov[1, 1]  # covYY(m2)
        PointsCov[point_id][5] = point_cov[1, 2]
        PointsCov[point_id][6] = point_cov[2, 0]
        PointsCov[point_id][7] = point_cov[2, 1]
        PointsCov[point_id][8] = point_cov[2, 2]  # covZZ(m2)
    if IdList:
        PointsCov_IdList = []
        for point_id in IdList:
            PointsCov_IdList.append(PointsCov[point_id])
        return PointsCov_IdList
    else:
        return PointsCov


def getPointsCov_analysis(chunk, IdList=[]):
    '''
    input:
        Metashape.document.chunk
        optional parameters:IdList is a set of point_id
            Note that the index of returned PointQuality are not align with point_id list but the given IdList.
    output:
        PointsCov = [PointCov,PointCov,...,PointCov]
            ,where PointCov = [stdP(mm), stdE(mm), stdM(mm), stdF(mm), stdE/stdF]
    This function must follow the step of 'optimization' with option 'estimating points covariance' checked.
    '''
    points = chunk.point_cloud.points
    M, T, R = MsTrans.getCoordTransMat(chunk)

    # stdP(mm), stdE(mm), stdM(mm), stdF(mm), stdE/stdF
    PointsCov = [[0, 0, 0, 0, 0] for i in range(len(points))]
    for point_id, point in enumerate(points):
        # skipping invalid points
        if not point.valid or point.cov:
            continue
        # Transform point covariance matrix from CLCS to PJCS
        point_cov = R * point.cov * R.t()
        # to get the max and min length of covariance components, we need to calculate eigenvalue by SVD method
        usv = point_cov.svd()  # Return type: tuple(Matrix,Vector,Matrix)
        ExtremeCOV = usv[1]
        PointsCov[point_id][0] = 1000 * math.sqrt(point_cov[0, 0] + point_cov[1, 1] + point_cov[2, 2])
        PointsCov[point_id][1] = math.sqrt(ExtremeCOV[0]) * 1000
        PointsCov[point_id][2] = math.sqrt(ExtremeCOV[1]) * 1000
        PointsCov[point_id][3] = math.sqrt(ExtremeCOV[2]) * 1000
        if math.sqrt(ExtremeCOV[2]) == 0:
            PointsCov[point_id][4] = 0
        else:
            PointsCov[point_id][4] = math.sqrt(ExtremeCOV[0]) / math.sqrt(ExtremeCOV[2])
    if IdList:
        PointsCov_IdList = []
        for point_id in IdList:
            PointsCov_IdList.append(PointsCov[point_id])
        return PointsCov_IdList
    else:
        return PointsCov


def getPointsData(chunk, Covariance=True, Quality=True, IdList=[]):
    '''
    get points coordinate, covariance, and quality as required
    output:
        PointsData = [PointData,PointData,...,PointData]
        PointData = [point_id, X(m), Y(m), Z(m), stdX(mm), stdY(mm), stdZ(mm), RE, RU, PA, IC]
    '''
    points = chunk.point_cloud.points
    crs = chunk.crs
    M, T, R = MsTrans.getCoordTransMat(chunk)

    # point_id, X(m), Y(m), Z(m), stdX(mm), stdY(mm), stdZ(mm), RE, RU, PA, IC
    PointsData = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for i in range(len(points))]
    for point_id, point in enumerate(points):
        # skipping invalid points
        if not point.valid:
            continue
        # Transform point coordinate from CLCS to PJCS
        point_coord = MsTrans.transPointCoord_CLCS2PJCS(point.coord, crs, M)
        # write result
        PointsData[point_id][0] = point_id
        PointsData[point_id][1] = point_coord[0]
        PointsData[point_id][2] = point_coord[1]
        PointsData[point_id][3] = point_coord[2]
    if Covariance:
        for point_id, point in enumerate(points):
            # skipping invalid points
            if not point.valid:
                continue
            # Transform point covariance matrix from CLCS to PJCS
            point_cov = R * point.cov * R.t()
            # write result
            PointsData[point_id][4] = 1000 * math.sqrt(point_cov[0, 0])
            PointsData[point_id][5] = 1000 * math.sqrt(point_cov[1, 1])
            PointsData[point_id][6] = 1000 * math.sqrt(point_cov[2, 2])
    if Quality:
        PointsCriterion = getPointsCriterion(chunk)
        for point_id, point in enumerate(points):
            # skipping invalid points
            if not point.valid:
                continue
            # Get point quality values
            PointCriterion = PointsCriterion[point_id]
            # write result
            PointsData[point_id][7] = PointCriterion[0]
            PointsData[point_id][8] = PointCriterion[1]
            PointsData[point_id][9] = PointCriterion[2]
            PointsData[point_id][10] = PointCriterion[3]
    if IdList:
        PointsData_IdList = []
        for point_id in IdList:
            PointsData_IdList.append(PointsData[point_id])
        return PointsData_IdList
    else:
        return PointsData
