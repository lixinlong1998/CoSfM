import numpy as np
from scipy.stats import gaussian_kde

'''
This script provide functions for analysing the spatial distribution of given points
Note that: the function need a third library: numpy, so the script use these function should be run in CMD
'''


def getBoundary(Points):
    '''
    input:
        Points={point_id:Point,...}
            ,where Point = [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, track_id]
                   Coordinates=[[X, Y, Z],[X, Y, Z]]   # in CLCS & PJCS respectively
                   Covariance = 3x3 covariance matrix
                   StandardDeviation = [stdX, stdY, stdZ]  # in mm
                   Quality = [RE, RU, PA, IC]
    output:
        Boundary_CLCS = [Points_Coord_Max[0], Points_Coord_Min[0], Points_Coord_Max[1], Points_Coord_Min[1]]
        Boundary_PJCS = [Points_Coord_Max[2], Points_Coord_Min[2], Points_Coord_Max[3], Points_Coord_Min[3]]
    '''
    Points_Coord = []
    for point_id, Point in Points.items():
        Coordinates_CLCS = Point[0][0]
        Coordinates_PJCS = Point[0][1]
        Points_Coord.append([Coordinates_CLCS[0], Coordinates_CLCS[1], Coordinates_PJCS[0], Coordinates_PJCS[1]])
    Points_Coord = np.asarray(Points_Coord)
    Points_Coord_Max = Points_Coord.max(axis=0)
    Points_Coord_Min = Points_Coord.min(axis=0)
    # boundary_X_MAX,boundary_X_MIN,boundary_Y_MAX,boundary_Y_MIN
    # right(east),left(west),upper(north),lower(south)
    Boundary_CLCS = [Points_Coord_Max[0], Points_Coord_Min[0], Points_Coord_Max[1], Points_Coord_Min[1]]
    Boundary_PJCS = [Points_Coord_Max[2], Points_Coord_Min[2], Points_Coord_Max[3], Points_Coord_Min[3]]
    return Boundary_CLCS, Boundary_PJCS


def assessDensity(PointList, crs='CLCS'):
    '''
    PointList = [Point,Point,...,Point]
        ,where Point = [Coordinates, Covariance, StandardDeviation, Quality, selected, valid, track_id]
               Coordinates=[[X, Y, Z],[X, Y, Z]]   # in CLCS & PJCS respectively
               Covariance = 3x3 covariance matrix
               StandardDeviation = [stdX, stdY, stdZ]  # in mm
               Quality = [RE, RU, PA, IC]
    '''
    Points_Coord = []
    if crs == 'CLCS':
        for Point in PointList:
            Coordinates_CLCS = Point[0][0]
            Points_Coord.append([Coordinates_CLCS[0], Coordinates_CLCS[1]])
    elif crs == 'PJCS':
        for Point in PointList:
            Coordinates_PJCS = Point[0][1]
            Points_Coord.append([Coordinates_PJCS[0], Coordinates_PJCS[1]])
    Points_Coord = np.asarray(Points_Coord)
    Points_Coord_X = Points_Coord[:, 0]
    Points_Coord_Y = Points_Coord[:, 1]
    values = np.vstack([Points_Coord_X, Points_Coord_Y])
    PointsDensity = gaussian_kde(values)(values)
    return PointsDensity
