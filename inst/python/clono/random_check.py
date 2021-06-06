import numpy as np
import matplotlib.pyplot as plt
import skimage.io
import pandas as pd
import re
import os
import glob

def remove_bg_plot (root, img_name, show_cluster=True, thres=0):
    '''
    On the left panel, plot the DAPI image in blue
    On the right panel, plot the final segmentation result super-imposed on the
    DAPI image
    Examples:
    remove_bg_plot (DNA_img, seg_large, seg_small)
    '''
    ori_img = skimage.io.imread (root+'/'+img_name)
    regex = re.compile ('_DAPI.tif')
    seg_img = skimage.io.imread (root+'/'+regex.sub ('_nuc_seg.tiff', img_name))

    if show_cluster:
        stat = pd.read_csv (root+'/'+regex.sub ('_nuc_seg.csv', img_name))
        cluster_n = stat[['cluster', 'condition']].groupby('cluster').aggregate({'condition':'count'})
        print (cluster_n)
        for i in range (len (stat)):
            if cluster_n.reset_index() ['condition'][stat['cluster'][i]+1] >= thres:
                cluster_lab = stat['cluster'][i]+1
            else: cluster_lab = 0
            seg_img [seg_img==stat['label'][i]] = cluster_lab

    lab_seg = skimage.color.label2rgb (seg_img, ori_img, bg_label=0, alpha=0.2)
    plt.subplot(1,2,1)
    zero_fill = np.zeros ([*ori_img.shape, 2])
    blue_img = np.concatenate ([zero_fill, ori_img [:,:,np.newaxis]], axis=-1)
    plt.imshow (blue_img/255)
    plt.axis ('off')

    plt.subplot(1,2,2)
    plt.imshow (lab_seg)
    plt.axis ('off')
    plt.show ()

root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/clono/clono/'
TF = 'nr2f2'
series='07'
remove_bg_plot (root+TF, '210317_siRNA_'+TF+'_00'+series+'_DAPI.tif', show_cluster=True)

def random_select_img (root_dir, save_dir, write_mode='a'):
    '''Randomly select images from the folders in the root directory
    The root directory contains a folder for each experimental condition.
    The naming convention of each image file is: $condition_$trial_DAPI.tif
    Examples:
    >>> random_select_img (root, root+'/proc_data', write_mode='w')
    '''
    all_files = glob.glob (root+'/*/*_DAPI.tif')
    rand_index = np.random.choice(len(all_files), 10)
    regex = re.compile ('^'+root_dir)
    regex2 = re.compile ('_DAPI.tif$')
    all_sel = [regex2.sub ('', regex.sub ('', all_files [i])) for i in rand_index]

    regex3 = re.compile ('/.*$')
    rand_df = pd.DataFrame.from_dict ({'condition': [regex3.sub ('', i) for i in all_sel]})
    regex4 = re.compile ('^.*_')
    rand_df ['trial'] =  [regex4.sub ('', i) for i in all_sel]
    rand_df ['true_num'] = 0

    save_name = save_dir+'/validation.csv'
    if write_mode=='a':
        header= False if os.path.exists (save_name) else True
    else: header = True
    rand_df.to_csv (save_name, mode=write_mode, header=header)

def clono_seg_plot (img_dir, img_file, stats_df, save_dir=None):
    '''
    Save 4 files in total:
        DAPI image in blue in tiff + pdf
        final segmentation result super-imposed on the DAPI image in tiff + pdf
    Args:
        stats_df: a pandas dataframe of the quantified (and filtered) statistics
    Examples:
    '''
    ori_img = skimage.io.imread (img_dir+'/'+img_file)
    if len (ori_img.shape) >2: ori_img = ori_img.sum (2)
    ori2seg = re.compile ('_DAPI.tif')
    seg_img_name = ori2seg.sub ('_nuc_seg.tiff', img_file)
    seg_img = skimage.io.imread (img_dir+'/'+seg_img_name)
    new_img = np.zeros (seg_img.shape)

    # select the appropriate rows corresponding to the image
    run_regex = re.compile ('run[0-9]+')
    Date_stat = run_regex.search (img_file).group(0)

    condition_regex = re.compile ('run[0-9]+_.*_')
    condition = condition_regex.search (img_file).group(0)
    condition = run_regex.sub ('', condition)
    condition_regex2 = re.compile ('^_')
    condition = condition_regex2.sub ('', condition)
    condition_regex3 = re.compile ('_.*$')
    condition = condition_regex3.sub ('', condition)

    series_regex = re.compile ('00[0-9]+_DAPI.tif$')
    series = series_regex.search (img_file).group(0)
    series_regex2 = re.compile ('_DAPI.tif$')
    series = int (series_regex2.sub ('', series))

    stats_fil = stats_df[(stats_df.Date.values == Date_stat) &
            (stats_df.condition.values ==condition.upper ()) & 
            (stats_df.series.values == series)]
    for i in stats_fil.index:
        new_img[seg_img==stats_fil.label[i]] = stats_fil.cluster[i] +1

    lab_seg = skimage.color.label2rgb (new_img, ori_img/255, bg_label=0)
    if save_dir is None: save_dir = img_dir
    save_file = ori2seg.sub ('_seg.png', img_file)
    plt.imsave (save_dir+save_file, lab_seg)

root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/clono/'
stats_df = pd.read_csv (root+'/combined_nuclei.csv')

fig_dir=root+'/figure_examples/GFP/'
save_dir=root+'/figure_save/GFP/'
for i in os.listdir (fig_dir): 
    if 'DAPI' in i: clono_seg_plot (fig_dir, i, stats_df, save_dir)

