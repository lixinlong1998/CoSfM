import math

'''
calculate statistic value of given list
'''


def removeZero(alist):
    output = []
    for i in alist:
        if i:
            output.append(i)
    return output


def removeRowsContainsZero(alist2d):
    output = []
    for row in alist2d:
        Count_zero = 0
        for i in row:
            if i == 0:
                Count_zero += 1
        if Count_zero == 0:
            output.append(row)
    return output


def removeRows(alist2d, usecols):
    '''
    usecols=[col,col,...]
    '''
    output = []
    for row in alist2d:
        Count_zero = 0
        for i in usecols:
            if row[i] == 0:
                Count_zero += 1
        if Count_zero == 0:
            output.append(row)
    return output


def calculateMAE(residuals):
    absResiduals = [abs(i) for i in residuals]
    return sum(absResiduals) / len(absResiduals)


def calculateRMSE(residuals):
    sqrResiduals = [i * i for i in residuals]
    return math.sqrt(sum(sqrResiduals) / len(sqrResiduals))

def calculateSTD(residuals):
    mean = sum(residuals) / len(residuals)
    sqrResiduals = [(i-mean)**2 for i in residuals]
    if len(sqrResiduals)>1:
        return math.sqrt(sum(sqrResiduals) / (len(sqrResiduals)-1))
    else:
        return math.sqrt(sum(sqrResiduals))

def calculateMedian(alist):
    alist.sort()
    half = len(alist) // 2
    return (alist[half] + alist[~half]) / 2


def calculatePercentile(alist, n):
    '''n = 1,2,3'''
    if n not in [1, 2, 3]:
        return False
    alist.sort()
    position = 1 + (len(alist) - 1) * n / 4
    pos_integer = int(math.modf(position)[1])
    pos_decimal = position - pos_integer
    quartile = alist[pos_integer - 1] + (alist[pos_integer] - alist[pos_integer - 1]) * pos_decimal
    return quartile


def listStatistic(alist):
    PointsCov_100 = max(alist)
    PointsCov_75 = calculatePercentile(alist, 3)
    PointsCov_50 = calculatePercentile(alist, 2)
    PointsCov_25 = calculatePercentile(alist, 1)
    PointsCov_0 = min(alist)
    PointsCov_avg = sum(alist) / len(alist)
    return [PointsCov_avg, PointsCov_0, PointsCov_25, PointsCov_50, PointsCov_75, PointsCov_100]
