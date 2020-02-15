# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 10:13:37 2020

@author: Mark
"""
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import join
import dicom
import pydicom
from pydicom.dataset import Dataset, FileDataset
from imageio import imread, imwrite
import astra

input_dir = r'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Data\11.2\DICOM 11.2\CT 11.2'
output_dir = r'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Data\11.2\DICOM 11.2\CT 11.2\recon'
dcm_file_list= ['IM-0001-0006.dcm','IM-0022-0027.dcm','IM-0023-0028.dcm','IM-0020-0025.dcm','IM-0019-0024.dcm']
dcm_file= dcm_file_list[0]
ds = pydicom.read_file(join(input_dir, dcm_file))
arr=ds.pixel_array

distance_source_origin = int(ds.DistanceSourceToEntrance)  # [mm]
distance_origin_detector = int(ds.DistanceSourceToDetector)  # [mm]
detector_pixel_size = float(ds.ImagerPixelSpacing[0])  # [mm]
detector_rows = int(ds.Rows)  # Vertical size of detector [pixels].
detector_cols = int(ds.Columns)  # Horizontal size of detector [pixels].
num_of_projections = int(ds.NumberOfFrames)
angles = np.linspace(0, 2 * np.pi, num=num_of_projections, endpoint=False)

projections = np.zeros((detector_rows, num_of_projections, detector_cols))
for i in range(num_of_projections):
    im = arr[i,:,:]
    im = im/65535
    projections[:, i, :] = im

# Copy projection images into ASTRA Toolbox.
proj_geom = \
  astra.create_proj_geom('cone', 1, 1, detector_rows, detector_cols, angles,
                         (distance_source_origin + distance_origin_detector) /
                         detector_pixel_size, 0)
projections_id = astra.data3d.create('-sino', proj_geom, projections)
 
# Create reconstruction.
vol_geom = astra.creators.create_vol_geom(detector_cols, detector_cols,
                                          detector_rows)

reconstruction_id = astra.data3d.create('-vol', vol_geom, data=0)
alg_cfg = astra.astra_dict('FP_CUDA')
alg_cfg['ProjectionDataId'] = projections_id
alg_cfg['ReconstructionDataId'] = reconstruction_id
#=============================================================================
algorithm_id = astra.algorithm.create(alg_cfg)
astra.algorithm.run(algorithm_id)
reconstruction = astra.data3d.get(reconstruction_id)
 
# Limit and scale reconstruction.
reconstruction[reconstruction < 0] = 0
reconstruction /= np.max(reconstruction)
reconstruction = np.round(reconstruction * 255).astype(np.uint8)
 
# Save reconstruction.
for i in range(detector_rows):
    im = reconstruction[i, :, :]
    im = np.flipud(im)
    imwrite(join(output_dir, 'reco%04d.png' % i), im)
 
# Cleanup.
astra.algorithm.delete(algorithm_id)
astra.data3d.delete(reconstruction_id)
astra.data3d.delete(projections_id)
