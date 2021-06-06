# This script determines the threshold of gene expression beyond which a cell
# is called positive

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skimage.io
import struct

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

root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/data/fkpd2/'
# for HLA-G
filter_img (root, '0', [5,10,15], gene='HLAG', clip_thres=50)
# for CGB
filter_img (root, '34', [6,7,8], gene='CGB', clip_thres=50)
# for TFAP2C
filter_img (root, '13', [4,5,7], gene='TFAP2C', clip_thres=50)

def show_IF (root, trial, series_num, vert_labels, channel_names=['HLAG', 'CGB',
    'TFAP2C', 'DAPI'], colors=['bf0504', '00c1ab', 'be9c00', '00acfc']):
    '''
    Show IF images
    Description: 
        Show all the channels per row
        Each row represents one condition
    Args:
        `root`: folder containing the image files
        `trial`: a list containing the experimental trials/batches
        `series_num`: a list containing the series to show
        `vert_labels`: experimental condition names
        `channel_names`: which channels to show
    '''
    # color coding
    color_prop = {one_chan:np.array (struct.unpack('BBB', 
        bytes.fromhex(one_c)))/255 for one_chan, one_c in zip (channel_names, colors)}

    fig, ax = plt.subplots (len (series_num), len (channel_names)+1)
    fig.set_size_inches ((len(channel_names)+1)*2, len(series_num)*2)
    for si, (one_trial, one_series, one_vert) in enumerate (zip(
        trial, series_num, vert_labels)):
        overlay_img = []
        for ci, one_chan in enumerate (channel_names):
            img_name = root+'/'+one_trial+'/series'+one_series+'_'+one_chan+'.tiff'
            img_bw = skimage.io.imread (img_name)
            if img_bw.shape[0] > img_bw.shape[2]:
                img_bw = np.transpose (img_bw, axes = [2,0,1])
            img_bw = np.squeeze (img_bw.max(0))
            proportion = color_prop [one_chan].reshape ([1,1,3])
            img_rgb = np.stack ([img_bw]*3, axis=-1)*proportion
            if one_chan != 'HLAG':
                norm_img = img_rgb/img_rgb.max()
            if one_chan == 'HLAG':
                norm_img = img_rgb.clip (0, 50)/50
            overlay_img.append (norm_img)
            ax[si, ci].imshow (norm_img)
            ax[si, ci].axis ('off')
            if si == 0:
                ax[si, ci].set_title (one_chan)
            if ci == 0:
                ax[si, ci].text(-300, 600, vert_labels[si])

        overlay_img = np.mean(np.stack (overlay_img, axis=-1), axis=-1)
        ax[si, ci+1].imshow (overlay_img/overlay_img.max())
        ax[si, ci+1].axis ('off')
        if si == 0: ax[si, ci+1].set_title ('overlay')

root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/data/'
save_dir='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/manuscript/'
show_IF (root, ['all', 'all', 'all'], ['33', '21', '13'], ['base', 'PD03', 'FK'])
plt.savefig (save_dir+'figure5/IF_example.pdf', dpi=700)
plt.show()
