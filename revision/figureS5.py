# validation for screening methods
import numpy as np
import matplotlib.pyplot as plt
import skimage.io
import pandas as pd
import re
import os
import glob

# ----------clonogenicity assay----------
def clono_seg_plot (stats_df, root, Date, condition, series, save_dir=None,
        Date_stat=None):
    '''
    Save 4 files in total:
        DAPI image in blue in tiff + pdf
        final segmentation result super-imposed on the DAPI image in tiff + pdf
    Args:
        stats_df: a pandas dataframe of the quantified (and filtered) statistics
    Examples:
    '''
    all_files = glob.glob (root+'/'+Date+'/*/*_nuc_seg.tiff')
    img_match = re.compile ('^.*/'+condition+'/.*'+series, re.IGNORECASE)
    seg_img_name = [i for i in all_files if img_match.match (i)]
    seg_img = skimage.io.imread (seg_img_name[0])
    new_img = np.zeros (seg_img.shape)

    seg2ori = re.compile ('_nuc_seg.tiff')
    ori_img_name = seg2ori.sub ('_DAPI.tif', seg_img_name[0])
    ori_img = skimage.io.imread (ori_img_name)

    # select the appropriate rows corresponding to the image
    if Date_stat is None: Date_stat = Date
    stats_fil = stats_df[(stats_df.Date.values==Date_stat) &
            (stats_df.condition.values ==condition.upper ()) & 
            (stats_df.series.values == int (series))]
    for i in stats_fil.index:
        new_img[seg_img==stats_fil.label[i]] = stats_fil.cluster[i] +1

    lab_seg = skimage.color.label2rgb (new_img, ori_img, bg_label=0, alpha=0.2)
    zero_fill = np.zeros ([*ori_img.shape, 2])
    blue_img = np.concatenate ([zero_fill, ori_img [:,:,np.newaxis]], axis=-1)

    if save_dir is not None:
        plt.imshow (blue_img/255)
        plt.axis ('off')
        plt.savefig (save_dir+'/clono_DAPI.pdf', dpi=700)
        plt.imsave (save_dir+'/clono_DAPI.png', blue_img/255)

        plt.imshow (lab_seg)
        plt.axis ('off')
        plt.savefig (save_dir+'/clono_seg.pdf', dpi=700)
        plt.imsave (save_dir+'/clono_seg.png', lab_seg)
    else: plt.imshow (lab_seg)
    plt.show()

root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/clono/'
save_dir='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/manuscript/figureS5/'
stats_df = pd.read_csv (root+'/filtered_nuclei.csv')
stats_df = pd.read_csv (root+'/combined_nuclei.csv')
#clono_seg_plot (stats_df, root, 'clono', 'ovol1', '08', save_dir, Date_stat='run1')
clono_seg_plot (stats_df, root, 'clono7', 'gfp', '12', None, Date_stat='run4')

# ----------signalling pathway screen----------
def cyto_seg_plot (root, series, save_dir):
    seg_org = skimage.io.imread (root+'/'+'series'+series+'_cyto_seg.tiff')
    cgb = skimage.io.imread (root+'/'+'series'+series+'_CGB.tiff')
    if cgb.shape[0] >cgb.shape[2]: cgb= np.transpose (cgb, [2,0,1])
    cgb = cgb.max (0)
    lab_seg = skimage.color.label2rgb (seg_org, cgb, bg_label=0, alpha=0.2)

    plt.imshow (cgb)
    plt.axis ('off')
    plt.savefig (save_dir+'/signal_CGB.pdf', dpi=700)
    plt.imsave (save_dir+'/signal_CGB.png', cgb)

    plt.imshow (lab_seg)
    plt.axis ('off')
    plt.savefig (save_dir+'/signal_cytoseg.pdf', dpi=700)
    plt.imsave (save_dir+'/signal_cytoseg.png', lab_seg)
    plt.show ()

root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/data/all/'
cyto_seg_plot (root, '13', save_dir)
