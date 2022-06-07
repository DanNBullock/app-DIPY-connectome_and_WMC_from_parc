#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 17:13:31 2022

@author: dan
"""


import os
import sys

sys.path.append('wma_pyTools')
startDir=os.getcwd()
#some how set a path to wma pyTools repo directory
#wmaToolsDir='wma_pyTools'
#wmaToolsDir='..'
#os.chdir(wmaToolsDir)
print(os.getcwd())
print(os.listdir())
import wmaPyTools.roiTools
import wmaPyTools.analysisTools
import wmaPyTools.segmentationTools
import wmaPyTools.streamlineTools
import wmaPyTools.visTools
import wmaPyTools.genUtils

from dipy.tracking.utils import reduce_labels
from dipy.tracking import utils

#os.chdir(startDir)
from scipy.io import savemat
import os
import json
import numpy as np
import nibabel as nib
import pandas as pd


# load inputs from config.json
with open('config.json') as config_json:
	config = json.load(config_json)


if not os.path.exists(os.path.join('wmc')):
    os.makedirs(os.path.join('wmc'))
if not os.path.exists(os.path.join('connectome')):
    os.makedirs(os.path.join('connectome'))
    
parcIn=nib.load(config['parc'])
lookupTable=wmaPyTools.genUtils.parcJSON_to_LUT(config['label'])

[renumberedAtlasNifti,reducedLookupTable]=wmaPyTools.analysisTools.reduceAtlasAndLookupTable(parcIn,lookupTable,removeAbsentLabels=True)

#load and oreint the streamlines
inTractogram=nib.streamlines.load(config['track'])
#try this, previously it may have been making inf values, but that was probably
#mrtrix misbehaving
orientedStreams=wmaPyTools.streamlineTools.orientAllStreamlines(inTractogram.streamlines)

M, grouping=utils.connectivity_matrix(orientedStreams, np.round(renumberedAtlasNifti.affine,2), label_volume=renumberedAtlasNifti.get_data().astype(int),
                        symmetric=False,
                        return_mapping=True,
                        mapping_as_streamlines=False)

#parse threshold
if 'threshold' in config.keys():
    if isinstance(config['threshold'],str):
        if config['threshold'].isnumeric():
            threshold=float(config['threshold'])
        else: 
            threshold=None
    elif isinstance(config['threshold'],int) or isinstance(config['threshold'],float):
        threshold=config['threshold']
    else:
        threshold=None   
else:
    threshold=None
    
#now determine if it is count or proportional
if not threshold== None:
    #if it's less than 1 its proportional
    if threshold<1:
        #we need to convert it to count
        thresholdCount=int(np.round(threshold*len(orientedStreams)))
        
    
    #otherwise it's count, do nothing
    else:
        thresholdCount=threshold
    
    #now iterate through 
    print ('Imposing streamline count threshold of ' +str(thresholdCount) )
    tempGrouping={}
    #probably a quicker and easier way to do this just using the M matrix and np.where
    for iConnections in list(grouping.keys()):
        #if there are a sufficient number of streams in it
        if len(grouping[iConnections])>=thresholdCount:
            #add it to the tempGrouping
                tempGrouping[iConnections]=grouping[iConnections]
        else:
            #do nothing, don't add it
            pass
    #set the output
    grouping=tempGrouping
    
    #now threshold the M, set the relevant entries to 0
    M[M<threshold]=0


#classification production
classification=wmaPyTools.streamlineTools.wmc_from_DIPY_connectome(grouping,reducedLookupTable)
#save it
savemat(os.path.join('wmc','classification.mat'),{ "classification": {"names": np.array(classification['names'], dtype=np.object), "index": classification['index'] }})
 
#create the conmat
wmaPyTools.genUtils.bl_conmat_fromDIPYandParc(M,reducedLookupTable,os.path.join('connectome'))

#create the network output
#conmat_to_JGFZ(arrayORcsv,indexIn,labelIn)
#no, just use the app, this format is rediculous
#https://brainlife.io/app/5f1fa182beafe9548f622a7e
#wmaPyTools.genUtils.conmat_to_JGFZ(M,os.path.join('connectome','index.json'),os.path.join('connectome','label.json'))
