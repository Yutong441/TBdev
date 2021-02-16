import numpy as np
import skimage.io
from skimage import measure
import matplotlib.pyplot as plt
import pandas as pd
from scipy import ndimage
from matplotlib.patches import Circle

def pixelcount(regionmask):
    return np.sum(regionmask)

def center_all_objects (mask):
    centroid = measure.regionprops (mask, extra_properties=(pixelcount,))
    N = len (centroid)
    centroid_arr = np.stack ([centroid[i].centroid for i in range (N)], axis=0)
    pixel_arr = np.stack ([centroid[i].pixelcount for i in range (N)], axis=0)
    return centroid_arr, pixel_arr

def check_centroid (segmented, labelled, level, circle_size, show_img=True,
        ax=None):
    labelled = np.round (labelled)
    if ax is None: fig,ax = plt.subplots()
    ax.imshow (segmented [:,:,level])
    for xx,yy,zz in zip(labelled[:,0], labelled[:,1], labelled[:,2]):
        if zz == level:
            circ = Circle((yy, xx), circle_size)
            ax.add_patch(circ)
    if show_img: plt.show ()


class Organ:
    def __init__ (self, inp_img, logging):
        '''
        `inp_img`: a labelled image.
        `logging`: a file object to store the history of filtering.
        '''
        self.lab_mask = inp_img
        self.logging = logging

    def get_centers (self):
        self.centroid, self.pcounts = center_all_objects (self.lab_mask)

    def get_index (self):
        self.index_id = np.unique (self.lab_mask)
        # 0 does not count
        self.index_id = self.index_id [self.index_id != 0]

    def filter_size (self, threshold):
        self.logging.write('filtering {} objects from {} \n'.format ( sum (self.pcounts <
            threshold ), len (self.pcounts)) )
        # filter away the objects with size < threshold
        self.get_index ()
        new_num = np.array ([self.index_id[i] for i in range (len (self.pcounts)) if 
            self.pcounts [i] < threshold])
        index = np.isin (self.lab_mask.reshape(-1), new_num)
        new_mask = index.reshape (self.lab_mask.shape) 
        print ('recalculating the centroids')
        self.lab_mask [new_mask] = 0
        self.centroid, self.pcounts = center_all_objects (self.lab_mask)
        self.get_index ()

    def show_centers (self, level=None, circle_size=10, max_level=9, num_col=3,
            disp=True):
        if level is None: level = list(np.arange (max_level) )
        if len (level) > max_level: level = level [:max_level]
        num_row = int (np.ceil (len (level)/num_col ))

        fig, axes = plt.subplots (num_col, num_row)
        ori_mask = self.lab_mask > 0
        for i, (one_lev, one_ax) in enumerate (zip (level, axes.ravel () )):
            one_ax.set_title ('level = '+str(one_lev))
            try:
                check_centroid (ori_mask, self.centroid, one_lev, circle_size,
                        show_img=False, ax=one_ax)
            except: pass
        if disp: plt.show ()

def stat_object (mask, properties = ('centroid', 'orientation',
    'major_axis_length', 'minor_axis_length', 'label', 'perimeter', 'area')):
    '''
    Compute the statistics for each segmented objects
    Args:
        `mask`: a labelled mask
        `properties`: what statistical property to calculate. They must be
        within the options of `skimage.measure.regionprops_table`
    '''
    all_df = []
    for i in range (mask.shape[2]):
        centroid = measure.regionprops_table (mask [:,:,i], properties= properties)
        centroid_df = pd.DataFrame (centroid)
        centroid_df['level'] = i
        all_df.append (centroid_df)
    return pd.concat (all_df, axis=0)

def label_2D (bw_mask):
    '''
    Label 3D images plane by plane.
    This is because the 3D labelling function may label disconnected elements.
    Args:
        `bw_mask`: a binary mask consisting of only 0 and 1
    '''
    max_index = 0
    plane_list=[]
    for i in range (bw_mask.shape[2]):
        one_plane = skimage.measure.label (bw_mask[:,:,i])
        new_plane = one_plane + max_index
        new_plane [one_plane == 0] = 0
        max_index += int (one_plane.max())
        plane_list.append (new_plane)
    return np.stack (plane_list, axis=-1)
