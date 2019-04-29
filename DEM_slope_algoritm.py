'''
Header
Name:       DEM_slope_algoritm.py
Authors:    Marc Meijer
Date:       19.02.2019
'''

#import
import math
import arcpy
import numpy as np
from copy import deepcopy
from pprint import pprint

#functions
def _completePath(workspace, subdir, nameList):
    for ix in range(len(nameList)):
        nameList[ix] = workspace + "/" + subdir + "/" + nameList[ix]
    return nameList

def _DEMtoNumPyArray(inGrid):
    # information on input raster
    in_arr = arcpy.RasterToNumPyArray(inGrid, nodata_to_value = -9999)

def _value2rank(arr):
    rank_array = np.zeros(arr.shape)
    sorted_unique = np.unique(arr)
    l = len(sorted_unique)
    for ix, v in enumerate(sorted_unique[::-1]):
        rank_array[np.where(arr==v)] = l-ix
    return rank_array

def _swapIndex(lst):
    lst[0], lst[1] = lst[1], lst[0]
    return lst

def _checkOrientation(arr, pstart, pend):
    change_vector = np.array(pstart) - np.array(pend)
    mode = -7999
    if (change_vector[0] > change_vector[1]) and change_vector[0] > 0:
        #moves bottom to top
        print 'moves bottom to top'
        mode = 0
        return mode
    elif (change_vector[0] < change_vector[1]) and (change_vector[0] < 0):
        #moves top to bottom
        print 'moves top to bottom'
        mode = 1
        return mode
    elif (change_vector[0] < change_vector[1]) and (change_vector[1] > 0):
        #moves right to left, transpose points and array
        print 'moves right to left'
        mode = 2
        return mode
    elif (change_vector[0] > change_vector[1]) and (change_vector[1] < 0):
        #moves left to right, transpose points and array
        print 'moves left to right'
        mode = 3
        return mode
        
def _makeFilter(pstart, arr):
    filter_arr = arr[pstart[0]-1:pstart[0]+1, pstart[1]-1:pstart[1]+2]
    return filter_arr

def _indexFilter(arr, pstart):
    pstart_arr = np.array(pstart)
    index_arr = np.array([
                        [[pstart_arr[0] -1, pstart_arr[1] - 1], [pstart_arr[0] -1, pstart_arr[1]],[pstart_arr[0] -1, pstart_arr[1] + 1]],
                        [[pstart_arr[0], pstart_arr[1] - 1], [pstart_arr[0], pstart_arr[1]],[pstart_arr[0], pstart_arr[1] + 1]]
                        ])
    return index_arr

def _addRim(arr):
    arr = np.pad(arr, (1,1), 'constant', constant_values = -99)
    return arr

def _adjustPoint(pnt):
    pnt = np.array(pnt) + 1
    #read this array into a list for right format
    pnt = pnt.tolist()
    return pnt
    
def _searchPoint(pstart, pend, arr):
    #defining the filter, this filter will look at all the neighbouring cells of the starting point
    #(neighbouring points have indeces NOT related to the overall array)
    arr = _addRim(arr)
    pstart = _adjustPoint(pstart)
    pend = _adjustPoint(pend)
    filter_arr = _makeFilter(pstart, arr)
    fstart = np.array([1,1])        #define index of starting point in filter

    #Calculate indices of the filter locations to be used in calculating the distances to the end point for every indices in filter
    filter_indices_arr = _indexFilter(filter_arr, pstart)

    #Calculate the distance (index-wise) between end point and points of the filter (returns np.array same size as filter)
    dist_to_pend_arr = np.linalg.norm(filter_indices_arr - pend, ord=2, axis=2)
    
    #Calculate distances from starting point to the neighbours
    dist_arr = np.array([[dist_d, dist_y, dist_d],[dist_x, 0, dist_x]])

    #Calculate the angles from starting point to the neighbours
    np.seterr(divide='ignore', invalid='ignore')    #to ignore NaN values
    angle_arr = np.arctan((filter_arr - filter_arr[1,1])/dist_arr) / ((1./180.)* math.pi)
    
    #Get index of point which has an angle lower than the user defined limit AND smallest distance to the end point 
    ind_plow = np.where((angle_arr < s_max) & (angle_arr > s_min))
    if len(ind_plow[0]) == 0:
        # too steep, no possibility
        return -9999
    elif len(ind_plow[0]) == 1:
        presult = np.append(ind_plow[0], ind_plow[1])
    elif len(ind_plow[0]) > 1:
        #more angles found, looking for shortest distance
        #Select from low angles the point with the shortest distance to end 
        rank_arr = _value2rank(dist_to_pend_arr)
        ind_pshort = np.where(rank_arr == min(rank_arr[ind_plow]))
        if len(ind_pshort[0]) == 1:
            presult = np.append(ind_pshort[0], ind_pshort[1])
        elif len(ind_pshort[0]) > 1:
            #distances to end point are equal
            #Select from possible angles the lowest angle
            ind_best = np.where(angle_arr == min(angle_arr[ind_pshort]))
            #both angles are equal as wel as the distances
            if len(ind_best[0]) > 1:
                 return -8999
            presult = np.append(ind_best[0], ind_best[1])
            
    #get the index of the new point in respect to the original array
    pnew_array = pstart - fstart + presult - 1

    #read this array into a list for right format
    pnew = pnew_array.tolist()
    return pnew

def _CalculateRouteNorth(pstart, pend, arr):
    pnts_dict = {tuple(pstart):arr[pstart[0]][pstart[1]]}
    
    while pstart[0]!= pend[0]:
        pstart = _searchPoint(pstart, pend, arr)

        #check new output if it contains the flags
        if pstart == -9999:
            print 'too steep no possibiliy'
            return pnts_dict

        elif pstart == -8999:
            print 'Possible angles are equally steep AND distance to end point are equal, can go both ways'
            return pnts_dict

        pnts_dict.update({tuple(pstart):arr[pstart[0]][pstart[1]]})
    return pnts_dict
    
def _CalculateRouteSouth(pstart, pend, arr):
    #invert start and end point
    pstart_inv = deepcopy(pend)
    pend_inv = deepcopy(pstart)
    pnts_dict = {tuple(pstart_inv):arr[pstart_inv[0]][pstart_inv[1]]}
    pstart = _searchPoint(pstart_inv, pend_inv, arr)
    pnts_dict.update({tuple(pstart):arr[pstart[0]][pstart[1]]})
    while pstart[0]!= pend_inv[0]:
        pstart = _searchPoint(pstart, pend_inv, arr)
        #check new output if it contains the flags
        if pstart == -9999:
            print 'too steep no possibiliy'
            return pnts_dict

        elif pstart == -8999:
            print 'Possible angles are equally steep AND distance to end point are equal, can go both ways'
            return pnts_dict
        pnts_dict.update({tuple(pstart):arr[pstart[0]][pstart[1]]})
    return pnts_dict

def _CalculateRouteWest(pstart, pend, arr):
    pnts_dict = {tuple(pstart):arr[pstart[0]][pstart[1]]}
    pend = _swapIndex(pend)
    while pstart[1]!= pend[0]:
        # turn arr and points to calculate route facing ,,north''
        pstart = _swapIndex(pstart)
        arr = arr.T
        pstart = _searchPoint(pstart, pend, arr)

        #check new output if it contains the flags
        if pstart == -9999:
            print 'too steep no possibiliy'
            return pnts_dict

        elif pstart == -8999:
            print 'Possible angles are equally steep AND distance to end point are equal, can go both ways'
            return pnts_dict

        #turn arr and points to get correct input facing ,,west''
        arr = arr.T
        pstart = _swapIndex(pstart)
        pnts_dict.update({tuple(pstart):arr[pstart[0]][pstart[1]]})
    return pnts_dict

def _CalculateRouteEast(pstart, pend, arr):
    #invert start and end point
    pstart_inv = deepcopy(pend)
    pend_inv = deepcopy(pstart)
    pnts_dict = {tuple(pstart_inv):arr[pstart_inv[0]][pstart_inv[1]]}
    pstart = _swapIndex(pstart_inv)
    
    while pstart[1]!= pend[1]:
        # turn arr and points to calculate route facing ,,north''
        arr = arr.T
        pstart = _searchPoint(pstart, pend, arr)

        #check new output if it contains the flags
        if pstart == -9999:
            print 'too steep no possibiliy'
            return pnts_dict

        elif pstart == -8999:
            print 'Possible angles are equally steep AND distance to end point are equal, can go both ways'
            return pnts_dict

        #turn arr and points to get correct input facing ,,west''
        arr = arr.T
        pstart = _swapIndex(pstart)
        pnts_dict.update({tuple(pstart):arr[pstart[0]][pstart[1]]})
        pstart = _swapIndex(pstart)
    return pnts_dict

'''
main
'''
#user defined values
dist_x = 10 
dist_y = 10
dist_d = math.sqrt(dist_x**2 + dist_y**2)
s_max = 8.0 
s_min = s_max * -1

arr = np.array([
    [5,5,5,6,7,9,15],
    [5,5,5,6,8,9,15],
    [5,5,6,9,9,13,15],
    [5,6,7,8,9,9,13],
    [5,6,7,8,6,8,13],
    [5,6,6,6,9,6,13],
    [5,7,7,8,10,9,15]], dtype=float)

'''
Testing all the possible cases
'''

#North
print arr
pstart = [6,3]
pend = [2,3]
mode = _checkOrientation(arr, pstart, pend)
if mode == 0:
    pts_route = _CalculateRouteNorth(pstart, pend, arr)
    pprint(pts_route)

#Sout
pstart = [2,3]
pend = [6,3]
mode = _checkOrientation(arr, pstart, pend)
if mode == 1:
    pts_route = _CalculateRouteSouth(pstart, pend, arr)
    pprint(pts_route)

#West
pstart = [3,6]
pend = [3,2]
mode = _checkOrientation(arr, pstart, pend)
if mode == 2:
    pts_route = _CalculateRouteWest(pstart, pend, arr)
    pprint(pts_route)

#East
pstart = [3,2]
pend = [3,6]
mode = _checkOrientation(arr, pstart, pend)
if mode == 3:
    pts_route = _CalculateRouteEast(pstart, pend, arr)
    pprint(pts_route)

#DEM
dist_x = 30 
dist_y = 30
dist_d = math.sqrt(dist_x**2 + dist_y**2)
s_max = 10.0 
s_min = s_max * -1
worksp = "E:\Marc_laptop_studie_ichthus\studie\Dreden\Python\Route_Final_Project"
inGridName = "DEM.tif"
inGrid = _completePath(worksp, "Grid", [inGridName])[0]

arr = arcpy.RasterToNumPyArray(inGrid, nodata_to_value = -9999)
pstart = [600,452]
pend = [200,460]
mode = _checkOrientation(arr, pstart, pend)
pts_route = _CalculateRouteNorth(pstart, pend, arr)
pprint(pts_route)


