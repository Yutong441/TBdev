import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skimage.io
from skimage import measure
import tifffile as tif
import cp_pipeline as cpp

# run pipeline
def segment_3d (img_3d):
    # obtain maximal intensity projection
    max_img = {'DNA': img_3d.max (0)}
    pipeline_filename = "segment_cyto.cppipe"
    workspace = cpp.run_pipeline(pipeline_filename, max_img)
    seg_nuc = workspace.object_set.get_objects ('IdentifyPrimaryObjects')
    seg_cyto = workspace.object_set.get_objects ('IdentifySecondaryObjects')
    # obtain image mask
    return seg_nuc.get_labels()[0][0], seg_cyto.get_labels()[0][0]

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('f', type=str, help='DAPI image folder')
    args = parser.parse_args()

    all_files = os.listdir(args.f)
    import re
    regex = re.compile ('series[0-9]+_DAPI.tiff')
    regex2 = re.compile ('_DAPI.tiff')
    for one_file in all_files:
        if regex.match (one_file):
            print ('segmenting '+one_file)
            DNA_img = skimage.io.imread (args.f+'/'+one_file)
            # if the channel dimension is the last
            if DNA_img.shape[2] < DNA_img.shape[1]:
                DNA_img = np.transpose (DNA_img, axes=[2, 0, 1])
            seg_nuc, seg_cyto = segment_3d (DNA_img)

            save_name1 = regex2.sub ('_nuc_seg.tiff', one_file)
            save_name2 = regex2.sub ('_cyto_seg.tiff', one_file)
            tif.imsave (args.f+'/'+save_name1, seg_nuc, compress=6)
            tif.imsave (args.f+'/'+save_name2, seg_cyto, compress=6)
