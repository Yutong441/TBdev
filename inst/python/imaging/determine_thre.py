# This script determines the threshold of gene expression beyond which a cell
# is called positive

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skimage.io
root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/data/all/'

def filter_img (root, series_num, threshold, gene='HLAG', struct='cyto',
        clip_thres=40, num_col=3):
    '''
    Args:
        `root`: root directory where all the images and quantification results
        are stored
        `series_num`: which series to read
        `threshold`: gene expression threshold, a list
        `gene`: fluorescence image of which gene
        `struct`: either 'cyto' or 'nuc' for cytoplasm and nuclear masks
        respectively
        `clip_thres`: clip the maximum fluoroscence value
        `num_col`: show the results in how many columns
    Returns:
        matplotlib showing the thresholded image
    '''
    seg_org = skimage.io.imread (root+'/'+'series'+series_num+'_'+struct+'_seg.tiff')
    hlag = skimage.io.imread (root+'/'+'series'+series_num+'_'+gene+'.tiff')
    if hlag.shape[2] < hlag.shape[1]: hlag = hlag.max(2)
    else: hlag = hlag.max(0)
    # read the quantification results
    fluoro = pd.read_csv (root+'/series'+series_num+'_.csv', index_col=[0])
    # keep those cells above a threshold

    for index, one_thres in enumerate (threshold):
        child_index = fluoro [fluoro[gene].values > one_thres]['child'].values
        mask = np.isin (seg_org, list (child_index) )
        lab_seg = skimage.color.label2rgb (seg_org*mask, hlag, bg_label=0)
        plt.subplot (np.ceil (len(threshold)/num_col), num_col, index+1)
        plt.imshow (lab_seg)
        plt.imshow (hlag.clip (0, clip_thres), alpha=0.2)
        plt.title ('{} threshold = {}'.format (gene, one_thres))
    plt.show()

# for HLA-G
filter_img (root, '17', [6,7,8], gene='HLAG', clip_thres=50)
# for CGB
filter_img (root, '40', [6,7,8], gene='CGB', clip_thres=50)
# for TFAP2C
filter_img (root, '34', [6,7,8], gene='TFAP2C', clip_thres=50)
