#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 17:13:31 2022

@author: dan
"""


import os
import sys

sys.path.append('wma_pyTools/')
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

import os
import json
import numpy as np
import nibabel as nib
import pandas as pd


# load inputs from config.json
with open('config.json') as config_json:
	config = json.load(config_json)

outDir='output'
if not os.path.exists(outDir):
    os.makedirs(outDir)
if not os.path.exists(os.path.join(outDir,'wmc')):
    os.makedirs(os.path.join(outDir,'wmc'))
if not os.path.exists(os.path.join(outDir,'connectome')):
    os.makedirs(os.path.join(outDir,'connectome'))
    
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

classification=wmaPyTools.streamlineTools.wmc_from_DIPY_connectome(grouping,lookupTable)

#what format is network output?

