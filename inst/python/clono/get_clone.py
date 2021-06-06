import numpy as np
import matplotlib.pyplot as plt
import skimage.io
from sklearn import cluster
from organoid import Organ, stat_object

def plot_cluster (stat, lab_img, ori_img):
    '''
    Examples:
    >>> plot_cluster (stats, seg_img, ori_img)
    '''
    plt.subplot (1,2,1)
    plt.imshow (ori_img)
    plt.subplot (1,2,2)
    for i in range (len (stat)):
        lab_img [lab_img==stat['label'][i]] = stat['cluster'][i]+1
    plt.imshow (lab_img)
    plt.show ()

if __name__ == '__main__':
    import argparse
    import os
    import re

    parser = argparse.ArgumentParser()
    parser.add_argument('f', type=str, help='DAPI image folder')
    args = parser.parse_args()

    all_files = os.listdir(args.f)
    regex = re.compile ('^.*_nuc_seg.tiff$')
    regex2 = re.compile ('_nuc_seg.tiff')
    log_f = open (args.f+'/'+'series_log.txt', 'w')

    for one_file in all_files:
        if regex.match (one_file):
            seg_img = skimage.io.imread (args.f+'/'+one_file)
            log_f.write ('image dimension= {}, {} \n'.format(*seg_img.shape))

            org = Organ (seg_img[:,:,np.newaxis], log_f)
            org.get_centers ()
            org.get_index ()
            stats = stat_object (org.lab_mask)

            clustering = cluster.DBSCAN (eps=30, min_samples=4).fit (
                    stats[['centroid-0', 'centroid-1']])
            stats['cluster'] = clustering.labels_
            stats['condition'] = args.f
            stats['series'] = one_file
            stats.to_csv (args.f+'/'+regex2.sub('_nuc_seg.csv', one_file))

    log_f.close()
