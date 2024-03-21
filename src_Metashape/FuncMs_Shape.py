import Metashape


def addShapes(chunk, shapesPathList):
    '''
    input:
        chunk
        shapesPathList = [shapePath,shapePath,...,shapePath] *.shp文件所在的文件夹中必须包含*.dbf和*.prj
    output:
        list of <class 'Metashape.Shape'>
    '''
    crs = chunk.crs
    # import new Shapes
    for shape_id, shapePath in enumerate(shapesPathList):
        chunk.importShapes(path=shapePath, replace=False, boundary_type=Metashape.Shape.NoBoundary,
                           format=Metashape.ShapesFormatSHP, columns='nxyzd', delimiter=', ', group_delimiters=False,
                           skip_rows=0, crs=crs)
        print(type(chunk.shapes[shape_id]))
        chunk.shapes[shape_id].label = 'Shape{}'.format(shape_id)
        print(chunk.shapes[shape_id].key)  # shape_id
    print(chunk.shapes.items())
    return chunk.shapes


def selectShape(shapes, shapesLabels):
    '''
    input:
        list of <class 'Metashape.Shape'>
        shapesLabels = [shapesLabel,shapesLabel,...,shapesLabel]
    output:
        shapes are selected.
    '''
    for shape in shapes:
        if shape.label in shapesLabels:
            shape.selected = 1
    return
