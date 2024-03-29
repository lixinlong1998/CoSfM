import os
import sys
import time
import Metashape
import csv

'''
导出检查点，返回的是一个tracks集合，每个检查点对应一个track，它包含了该点在对应像片子集上的图像坐标。

CheckPoints = [CheckPoint,CheckPoint,...,CheckPoint]
CheckPoint = [marker_id, track]
track = [identifier, [projection.coord[0], projection.coord[1]]]
'''


def getCameraIdentifier(camera):
    '''
    create identifier
    '''
    # '2021:04:04 13:40:13' -> 20210404134013
    identifier = ''
    for s in camera.photo.meta['Exif/DateTimeOriginal']:
        if s.isdigit():
            identifier += s
    # 'DSC00009.JPG' -> 09  (# camera.label is DSC00004 in Metashape 1.8.5)
    # identifier += camera.label[-6:-4]
    identifier += camera.label[-2:]
    return int(identifier)

# -------------------------------------------------------------------------------- Loading Metashape document and chunks
starttime = time.perf_counter()

# access document
doc = Metashape.app.document

# choose the chunk enabled and with maximum images
chunk_images_num = 0
for chunk_i in Metashape.app.document.chunks:
    if chunk_i.enabled:
        chunki_images_num = len(chunk_i.cameras)
        if chunki_images_num >= chunk_images_num:
            chunk_images_num = chunki_images_num
            chunk = chunk_i
        else:
            continue

# creat a file for writing results
path_CPsdatabase = doc.path[0:(len(doc.path) - 4)] + '_CPsdatabase.txt'
reportFile = open(path_CPsdatabase, "w")
fwriter = csv.writer(reportFile, delimiter='\t', lineterminator='\n')

for marker_id, marker in enumerate(chunk.markers):
    print('[Script]    Marker:', marker_id)
    # skip ground control points(GCPs)
    if marker.reference.enabled:
        print('[Script]    Skip:', marker_id)
        continue
    track = []
    for camera, projection in marker.projections.items():
        identifier = getCameraIdentifier(camera)
        if projection.valid:
            if projection.pinned:
                track.append([str(identifier), [projection.coord[0], projection.coord[1]]])  # identifier is string
            else:
                print('blueflag')
        else:
            print('greyflag')
    print('[Script]    image count:', len(track))
    fwriter.writerow([marker_id, marker.label, track])
reportFile.close()

print('[Script][TimeCost]    Calculate CPs:', time.perf_counter() - starttime)

starttime = time.perf_counter()
print('[Script][TimeCost]    :', time.perf_counter() - starttime)