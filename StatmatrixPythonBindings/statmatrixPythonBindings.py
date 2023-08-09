import os
import platform
import time
from datetime import datetime
from tqdm import tqdm

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.io import loadmat
import mat73

import pickle

# -----------------------------------------------------------

""" 
Code to use the matlab stat matrix and behavior matrix files in python 
See also: https://github.com/FortinLab/statMatrixNeuralAnalysis

Overly verbose with the intention of toolbox use and easy adaptation of submodules

Code was writen with maintence in mind (and a lot of checks) but stay up to date with any updates

STAT MATRIX CONTENTS

    'statMatrix', 
    'statMatrixColIDs', 
    'unitSummary'

STAT MATRIX COLUMN IDs (view_statMat_cols)

    0 TimeBin
    1 T04_LFP_Raw
    2 T04_LFP_Raw_HilbVals
    3 T04_LFP_Theta
    4 T04_LFP_Theta_HilbVals
    5 T04_LFP_LowBeta
    6 T04_LFP_LowBeta_HilbVals
    7 T04_LFP_Beta
    8 T04_LFP_Beta_HilbVals
    9 T04_LFP_LowGamma
    10 T04_LFP_LowGamma_HilbVals
    11 T04_LFP_HighGamma
    12 T04_LFP_HighGamma_HilbVals
    13 T04_LFP_Ripple
    14 T04_LFP_Ripple_HilbVals
    15+ <Unit names>


Keiland Cooper Last updated, 210820
Fortin lab at UCI
"""


def path_constructor(dataPath, fnames): 
    """ """
    return [os.path.join(dataPath, fn) for fn in fnames]

def dataloader_behavior_matrix(matPath):
    """ """
    return loadmat(matPath)

def dataloader_stat_matrix_tetrodewise(matPath): 
    """ """
    return mat73.loadmat(matPath)

def parse_behavior_matrix(behavPaths, computeMeta=None):
    """  """
    behavMat = dataloader_behavior_matrix(behavPaths[0])
    
    behaviorMatrix = {}
    behaviorMatrix['taskMetadata'] = behavMat['behavMatrix']
    behaviorMatrix['colIds'] = behavMat['behavMatrixColIDs']
    
    behaviorMatrix['datetimeCompute'] = datetime.now().strftime("%y%m%d %H:%M:%S")
    behaviorMatrix['locationCompute'] = platform.node()
    behaviorMatrix['Path'] = os.environ['PWD']
    
    if computeMeta != None:
        for k in computeMeta.keys():
            behaviorMatrix[k] = computeMeta[k]
            
    return behaviorMatrix


def statmatrix_checks(tetData, stat_mat): 
    """ Run checks to make sure the code is not depriciated """
    # check to make sure the shapes are correct
    nRunits = tetData['unitData'].shape[1]
    if nRunits != tetData['nUnits']:
        raise ValueError(f"Reported number of units ({tetData['nUnits']}) does " \
                           f"not equal the captured number of units spike data ({nRunits})")

    # check to make sure the statmatrix cols haven't changed
    if (tetData['colLen']-tetData['nUnits']) != statMatColIdx['nStandardCols']:
        raise ValueError('Reported number of statmatrix cols does not equal ' \
                         'computed number. Please check & update "statMatColIdx"')

    # check to make sure the statmatrix col unit names are the same as the computed ones
    if stat_mat['statMatrixColIDs'][statMatColIdx['nStandardCols']:] != tetData['unitNames']:
        raise ValueError('Statmatrix col unit names are not the same as the unitSummary names')

# Store the statmat col idx in a 
# dict for easy maintenance
statMatColIdx = {'TimeBin':0,
                 'LFP_Raw':1,
                 'nStandardCols':15}

def compile_tetrode_data(stat_mat, statMatColIdx):
    """ Parse a single tetrode's statmatrix file """
    
    tetData = {}
    tetData['tetName'] = stat_mat['statMatrixColIDs'][1]

    # Collect the Time Data
    tetData['timeBin'] = stat_mat['statMatrix'][:,statMatColIdx['TimeBin']]

    # Collect the LFP Data
    tetData['lfpRaw'] = stat_mat['statMatrix'][:,statMatColIdx['LFP_Raw']]

    # Collect the Unit Data
    tetData['unitNames'] = stat_mat['statMatrixColIDs'][statMatColIdx['nStandardCols']:] #stat_mat['unitSummary']['UnitName']
    tetData['nUnits'] = len(tetData['unitNames'])
    tetData['unitData'] = stat_mat['statMatrix'][:,statMatColIdx['nStandardCols']:]
    tetData['colLen'] = len(stat_mat['statMatrixColIDs'])

    tetData['TimeComputed'] = datetime.now().strftime("%y%m%d %H:%M:%S")
    
    return tetData
    
def print_tetrode_info(tetData):
    print(f"Tetrode: {tetData['tetName']}")
    print(f"nUnits: {tetData['nUnits']}")
    print(f"{tetData['unitNames']}")
    print(f"{tetData['TimeComputed']}")
    

def build_multi_tetrode_objects(tetrodePaths, statMatColIdx=None, verbose_tetInfo=False, 
                                verbose_progressbar=False):
    """ 
    Builds tetrode objects from a list of given file names 
    A bit slow, takes a couple minutes to run on large datasets
    
    Params:
        tetrodePaths : list
            list of filenames of tetrodes to parse
        statMatColIdx : dict
            optional dictionary to update indicies
        verbose_tetInfo : bool
            optional print tetrode info
        verbose_progressbar : bool
            optional progress bar
    
    Returns:
        tetDataAll : list
            List of tet objects
    
    """
    
    # Store the statmat col idx in a 
    # dict for easy maintenance
    if statMatColIdx == None:
        statMatColIdx = {'TimeBin':0,
                         'LFP_Raw':1,
                         'nStandardCols':15}
    
    nTetrodes = len(tetrodePaths)
    tetDataAll = []
    for tetInd in tqdm(range(nTetrodes), disable=not verbose_progressbar):
        # Load in the data
        stat_mat = dataloader_stat_matrix_tetrodewise(tetrodePaths[tetInd])
        
        # Compile tet objects
        tetData = compile_tetrode_data(stat_mat, statMatColIdx)
        if verbose_tetInfo: print_tetrode_info(tetData)

        # Run checks
        statmatrix_checks(tetData, stat_mat)

        tetDataAll.append(tetData)
    
    return tetDataAll

def view_statMat_cols(stat_mat):
    """ """
    i = 0
    for c in stat_mat['statMatrixColIDs']:
        print(i, c)
        i += 1

def compute_total_units(tetDataAll):
    """ Compute total units across tetrodes """

    totalUnits = 0
    for tet in tetDataAll:
        totalUnits += tet['nUnits'] 
    
    return totalUnits

def build_lfp_matrix(tetDataAll):
    """ """
    lfps = []
    for tet in tetDataAll:
        lfps.append(tet['lfpRaw'])
    return np.swapaxes(np.array(lfps),0,1)

def build_unit_matrix(tetDataAll):
    """ """
    units = []
    for tet in tetDataAll:
        if tet['unitData'].shape[1] > 0: # drop empty units
            units.append(tet['unitData'])

    return np.concatenate(units, axis=1)

def build_tet_unit_names(tetDataAll):
    """ """
    tetNames, unitNames = [],[]
    for tet in tetDataAll:
        tetNames.append(tet['tetName'])
        unitNames.append(tet['unitNames'])

    return np.array(tetNames), np.concatenate(unitNames)

def parse_stat_matrix(tetrodePaths, computeMeta=None, statMatColIdx=None, 
                      verbose_tetInfo=False, verbose_progressbar=False):
    """ Compiles pythonic arrays of single file tetrode stat matrix files
        Wraps build_multi_tetrode_objects()
    
    Params:
        tetrodePaths : list
            list of filenames of tetrodes to parse
        statMatColIdx : dict
            optional dictionary to update indicies
        verbose_tetInfo : bool
            optional print tetrode info
        verbose_progressbar : bool
            optional progress bar
    
    Returns:
        statMatrix : dict
            compiled stat matrix files in numpy arrays with keys:
                lfpMat = matrix of continuous lfp data (time x channel)
                unitMat = matrix of continuous unit data (time x unit)
                tetNames = Names of tetrodes
                unitNames = Names of units
                totalUnits = Total number of units
                totalElectrodes = Total number of electrodes

                datetimeCompute = %y%m%d %H:%M:%S
                locationCompute = computer name
                computeDuration = computation duration
                Path = working directory
        
                    
    """
    
    timer_start = time.time()
    
    # Build tet objects
    tetDataAll = build_multi_tetrode_objects(tetrodePaths, statMatColIdx=statMatColIdx, 
                            verbose_tetInfo=verbose_tetInfo, 
                            verbose_progressbar=verbose_progressbar)
    
    
    # Format the data from the objects 
    statMatrix = {}
    statMatrix['lfpMat'] = build_lfp_matrix(tetDataAll)
    statMatrix['unitMat'] = build_unit_matrix(tetDataAll)
    statMatrix['tetNames'], statMatrix['unitNames'] = build_tet_unit_names(tetDataAll)
    statMatrix['totalUnits'] = compute_total_units(tetDataAll)
    statMatrix['totalElectrodes'] = len(tetDataAll)
    
    statMatrix['datetimeCompute'] = datetime.now().strftime("%y%m%d %H:%M:%S")
    statMatrix['locationCompute'] = platform.node()
    statMatrix['computeDuration'] = (time.time() - timer_start) / 60
    statMatrix['Path'] = os.environ['PWD']
    
    if computeMeta != None:
        for k in computeMeta.keys():
            statMatrix[k] = computeMeta[k]
            
    return statMatrix
    
def save_obj(data, file, binary=False):
    """ Save an object to a pickle file
    
        pickle.HIGHEST_PROTOCOL is binary file, else text file
    """
    if binary: 
        with open(file+'.pkl', 'wb') as f:
            pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
    else: 
        with open(file+'.pkl', 'wb') as f:
            pickle.dump(data, f)

def load_obj(file):
    """ """
    with open(file+'.pkl', 'rb') as f:
        return pickle.load(f)





if __name__ == "__main__":
	# Example code 
	
	# Compile data and path locations. Can be integrated elsewhere
	dataPath = os.path.join('..', 'CA3_data')

	behaviorfNames = ["HC01_071217_task-GEcut03_BehaviorMatrix.mat"]

	tetrodefNames = ["HC01_071217_task-GEcut03_T01.mat",
	                "HC01_071217_task-GEcut03_T02.mat",
	                "HC01_071217_task-GEcut03_T03.mat",
	                "HC01_071217_task-GEcut03_T04.mat",
	                "HC01_071217_task-GEcut03_T05.mat",
	                "HC01_071217_task-GEcut03_T06.mat",
	                "HC01_071217_task-GEcut03_T07.mat",
	                "HC01_071217_task-GEcut03_T08.mat",
	                "HC01_071217_task-GEcut03_T09.mat",
	                "HC01_071217_task-GEcut03_T11.mat",
	                "HC01_071217_task-GEcut03_T12.mat",
	                "HC01_071217_task-GEcut03_T13.mat",
	                "HC01_071217_task-GEcut03_T14.mat",
	                "HC01_071217_task-GEcut03_T15.mat",
	                "HC01_071217_task-GEcut03_T16.mat",
	                "HC01_071217_task-GEcut03_T17.mat",
	                "HC01_071217_task-GEcut03_T18.mat",
	                "HC01_071217_task-GEcut03_T20.mat",
	                "HC01_071217_task-GEcut03_T21.mat",
	                "HC01_071217_task-GEcut03_T22.mat",
	                "HC01_071217_task-GEcut03_T24.mat"]

	# Add useful metadata to the computed data structures
	computeMeta = {'subjName':'HC01', 'sessionName':'071217_task-GEcut03', 'authorName':'Keiland'}

	# Build the paths used by the functions and load in the data
	behaviorPaths = path_constructor(dataPath, behaviorfNames)
	behaviorMatrix = parse_behavior_matrix(behaviorPaths, computeMeta=computeMeta)

	tetrodePaths = path_constructor(dataPath, tetrodefNames)
	statMatrix = parse_stat_matrix(tetrodePaths, computeMeta=computeMeta, statMatColIdx=None, 
	                      verbose_tetInfo=False, verbose_progressbar=True)

	# Save the data to a file
	savePath = os.path.join(dataPath, 'statMatrix')
	save_obj(statMatrix, savePath)

