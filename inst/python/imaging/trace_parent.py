import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def find_parent (child, parent, index, threshold=0.1):
    child_index = child == index
    item_vec = parent [child_index]
    arr_ind, arr_count = np.unique (item_vec, return_counts=True)
    max_ind = arr_ind [arr_count == arr_count.max() ]
    if arr_count.max ()/child_index.sum() > threshold:
        return max_ind
    else:
        print ('selected item does not have a parent')
        return np.NaN

def find_parent_for_all (child, parent):
    '''
    Find which each object in `parent` that each object in `child` belongs to
    Args:
        `child`: a 3D image of smaller objects
        `parent`: a 3D image of larger objects
    '''
    all_ids = np.unique (child).astype (int)
    N = len (all_ids)
    obj_tree = np.zeros ([N, 2])
    obj_tree [:,0] = all_ids
    for i in range (N):
        obj_tree [i, 1] = int (find_parent (child, parent, obj_tree[i,0])[0])
    return obj_tree

def find_parent_fast (child, parent):
    child_flat = child.lab_mask.reshape (-1)
    paren_flat = parent.lab_mask.reshape (-1)
    ind = child_flat !=0
    apixel = np.stack ([child_flat[ind], paren_flat[ind]], axis=1)
    apixel = pd.DataFrame (apixel, columns=['child', 'parent'])
    apixel = apixel.groupby('child').agg(lambda x:x.value_counts().index[0])
    return apixel.reset_index()

def check_parent (mask_img, obj_tree, parent_img=None):
    new_mask = mask_img
    for i in range (len (obj_tree)):
        new_mask [new_mask == obj_tree [i,0] ] = obj_tree [i,1]
    plt.imshow (new_mask)
    if parent_img is not None: plt.imshow (parent_img, alpha=0.3, cmap='gray')
    plt.show ()

def get_fluoro (org_ob, fluoro, tree_arr, label):
    '''
    Args:
        `org_ob`: an Org object
        `fluoro`: fluorescence image to quantify
        `tree_arr`: a pandas dataframe
    '''
    N = int (org_ob.lab_mask.max())
    tree_arr [label] =0.
    for index, i in enumerate (tree_arr.child.astype(int)):
        # no need to minus or subtract indices because the `lab_mask` numbering
        # also starts from 1
        bw_mask = org_ob.lab_mask==i
        tree_arr [label].iloc[index] = np.sum (fluoro*bw_mask)/np.sum(bw_mask)
    return tree_arr

def get_fluoro_fast (org_ob, fluoro, tree_arr, label):
    # subtract away the background
    fluoro = (fluoro - fluoro.mean ()).clip (0, np.infty)
    org_flat = org_ob.lab_mask.reshape (-1)
    fluo_flat = fluoro.reshape (-1)
    ind = org_flat !=0
    org_fluo = np.stack ([org_flat [ind], fluo_flat[ind]], axis=1)
    org_fluo = pd.DataFrame (org_fluo , columns=['ID', 'FL'])
    mean_fl = org_fluo.groupby ('ID').mean ()
    mean_fl['ID'] = mean_fl.index
    tree_arr [label] = np.array ([mean_fl.FL.values[mean_fl.ID==i] for i in
            tree_arr.child.astype(int)])
    return tree_arr
