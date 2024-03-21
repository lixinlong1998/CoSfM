def extractDenseCloud(chunk, shapes):
    '''
    input:
        Metashape.document.chunk
        Metashape.Shape
    output:
        PointsCriterion = [PointCriterion,PointCriterion,...,PointCriterion]
            ,where PointCriterion = [RE, RU, PA, IC]
    Select dense points based on shapes.
    '''

    points = chunk.model.points.selectPointsByShapes([shape])

    return PointsCriterion
