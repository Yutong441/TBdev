import cp_pipeline as cpp

def segment_2d (img_2d):
    '''
    Segment the DAPI channel.
    Args:
        `img_2d`: a numpy array in the form of (H, W), which stands for
        height and width respectively. This can be an image of the DAPI channel
    Returns:
        a list containing 
        1. an array of shape (H, W) as the background segmentation
        1. an array of shape (H, W) as the background + nuclei segmentation 
    '''
    # obtain maximal intensity projection
    max_img = {'DNA': img_2d}
    pipeline_filename = "segment_background.cppipe"
    workspace = cpp.run_pipeline(pipeline_filename, max_img)
    seg_large = workspace.object_set.get_objects ('large_mask')
    seg_small = workspace.object_set.get_objects ('small_mask')
    # obtain image mask
    return seg_large.get_labels()[0][0], seg_small.get_labels()[0][0]

if __name__ == '__main__':
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import skimage.io
    from skimage import measure
    import tifffile as tif
    import argparse
    import re

    parser = argparse.ArgumentParser()
    parser.add_argument('f', type=str, help='DAPI image folder')
    args = parser.parse_args()

    all_files = os.listdir(args.f)
    regex = re.compile ('^.*_DAPI.tif$')
    regex2 = re.compile ('_DAPI.tif')

    for one_file in all_files:
        if regex.match (one_file):
            print ('segmenting '+one_file)
            DNA_img = skimage.io.imread (args.f+'/'+one_file)
            # if the image is 3D
            if len(DNA_img.shape) >2: DNA_img  = DNA_img.sum(2)
            seg_large, seg_small = segment_2d (DNA_img)

            save_name = regex2.sub ('_nuc_seg.tiff', one_file)
            large_bg = seg_large==0
            tif.imsave (args.f+'/'+save_name, seg_small*large_bg, compress=6)
