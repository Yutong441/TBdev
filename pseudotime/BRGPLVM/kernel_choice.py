'''
This script documents the design of kernels for different models.
Note that in gpflow1, when multiple kernels are defined in a function and
returned as a list, sometimes this function does not work.
Therefore, I have to directly define the kernels in the main script.

Furthermore, if possible, in gpflow1, use `gpflow.ekernels` instead of
`gpflow.kernels` . This is because in `gpflow.ekernels`, the psi statistics have
already been coded in a closed form, which makes computation faster. However,
this is not true for `gpflow.kernels`, for which psi statistics are evaluated
by quadrature.

This is again not an issue with gpflow2, which automatically detects
implementation of analytical solution. However, my feeling is that for simple
things such as zero kernel, quadrature is still used, which is completely
unnecessary. There is no way to change that. This is why I prefer gpflow1.
'''
import gpflow 

# ----------change point kernel with tie----------
# main branch models development from cMor to ICM, TB and CTB
# 1st branch (kernel1) models EVT
# 2nd branch (kernel2) models STB main
# 3rd and 4th branches model 2 subpopulations of STB
# reassign labels to reflect the hierarchy
branch_num = np.unique (branch.values)
branch = branch.values.reshape ([-1, 1]) == branch_num.reshape ([1, -1])
branch_inp = branch
branch_inp [:,0] = 1
branch_inp [:,2] = branch_inp [:,2] + branch_inp [:,3] +branch_inp [:,4]

# define the RBF components
kernel0 = gpflow.ekernels.RBF (1, variance=0.1)
kernel1 = gpflow.ekernels.RBF (1, variance=0.1)
kernel2 = gpflow.ekernels.RBF (1, variance=0.1)
kernel3 = gpflow.ekernels.RBF (1, variance=0.1)
kernel4 = gpflow.ekernels.RBF (1, variance=0.1)

# In order to constrain the branching point, it is important to specify kernel
# parameters before defining the kernel. Note this is not a problem with
# gpflow2 because of eager execution
BP1 = gpflow.params.Parameter ([1.])
BP2 = gpflow.params.Parameter ([2.])

cp_kernel1 = CP.tanh_CP ( [bu.zero_kern (), kernel1], locations=BP1)
cp_kernel2 = CP.tanh_CP ( [bu.zero_kern (), kernel2], locations=BP1)
cp_kernel3 = CP.tanh_CP ( [bu.zero_kern (), kernel3], locations=BP2)
cp_kernel4 = CP.tanh_CP ( [bu.zero_kern (), kernel4], locations=BP2)

kern_list = gpflow.ParamList([kernel0, cp_kernel1, cp_kernel2, cp_kernel3, cp_kernel4])

# ----------change point kernel no ties----------
# define the RBF components
kernel0 = gpflow.ekernels.RBF (1, variance=0.1)
kernel1 = gpflow.ekernels.RBF (1, variance=0.1)
kernel2 = gpflow.ekernels.RBF (1, variance=0.1)
kernel3 = gpflow.ekernels.RBF (1, variance=0.1)
kernel4 = gpflow.ekernels.RBF (1, variance=0.1)

BP1 = gpflow.params.Parameter ([1.])
BP2 = gpflow.params.Parameter ([1.])
BP3 = gpflow.params.Parameter ([2.])
BP4 = gpflow.params.Parameter ([2.])

cp_kernel1 = CP.tanh_CP ( [bu.zero_kern (), kernel1], locations=BP1)
cp_kernel2 = CP.tanh_CP ( [bu.zero_kern (), kernel2], locations=BP2)
cp_kernel3 = CP.tanh_CP ( [bu.zero_kern (), kernel3], locations=BP3)
cp_kernel4 = CP.tanh_CP ( [bu.zero_kern (), kernel4], locations=BP4)

kern_list = gpflow.ParamList([kernel0, cp_kernel1, cp_kernel2, cp_kernel3, cp_kernel4])

# ----------Change Point with constraint sequence of branching points----------
kernel0 = gpflow.ekernels.RBF (1, variance=0.1)
kernel1 = gpflow.ekernels.RBF (1, variance=0.1)
kernel2 = gpflow.ekernels.RBF (1, variance=0.1)
kernel3 = gpflow.ekernels.RBF (1, variance=0.1)
kernel4 = gpflow.ekernels.RBF (1, variance=0.1)

BP1 = gpflow.params.Parameter ([1.])
# constrain the BP2 to be later than BP1
diff = gpflow.params.Parameter ([1.], transform=gpflow.transforms.positive)

cp_kernel1 = CP.tanh_CP ( [bu.zero_kern (), kernel1], locations=BP1)
cp_kernel2 = CP.tanh_CP ( [bu.zero_kern (), kernel2], locations=BP1)
cp_kernel3 = CP.tanh_CP ( [bu.zero_kern (), kernel3], locations=BP1, constraints=diff)
cp_kernel4 = CP.tanh_CP ( [bu.zero_kern (), kernel4], locations=BP1, constraints=diff)

kern_list = gpflow.ParamList([kernel0, cp_kernel1, cp_kernel2, cp_kernel3, cp_kernel4])

# ----------RBF only----------
branch_num = np.unique (branch.values)
branch = branch.values.reshape ([-1, 1]) == branch_num.reshape ([1, -1])
branch_inp = branch
branch_inp [:,0] = 1

kernel0 = gpflow.ekernels.RBF (1, variance=0.1)
kernel1 = gpflow.ekernels.RBF (1, variance=0.1)
kernel2 = gpflow.ekernels.RBF (1, variance=0.1)
kernel3 = gpflow.ekernels.RBF (1, variance=0.1)
kernel4 = gpflow.ekernels.RBF (1, variance=0.1)
kern_list = gpflow.ParamList([kernel0, kernel1, kernel2, kernel3, kernel4])

# ----------3 branches----------
# for 'model_dict_3b'
branch_num = np.unique (branch.values)
branch = branch.values.reshape ([-1, 1]) == branch_num.reshape ([1, -1])
branch_inp = branch
branch_inp [:,0] = 1
branch_inp [:,2] = branch_inp [:,2] + branch_inp [:,3] +branch_inp [:,4]
branch_inp =branch_inp [:,:3] 

kernel0 = gpflow.ekernels.RBF (1, variance=0.1)
kernel1 = gpflow.ekernels.RBF (1, variance=0.1)
kernel2 = gpflow.ekernels.RBF (1, variance=0.1)
BP1 = gpflow.params.Parameter ([1.])
cp_kernel1 = CP.tanh_CP ( [bu.zero_kern (), kernel1], locations=BP1)
cp_kernel2 = CP.tanh_CP ( [bu.zero_kern (), kernel2], locations=BP1)
kern_list = gpflow.ParamList([kernel0, cp_kernel1, cp_kernel2])

# ----------Equal assignment to main branch----------
# for 'model_dict_3b_prob'
branch_num = np.unique (branch.values)
branch = branch.values.reshape ([-1, 1]) == branch_num.reshape ([1, -1])
branch_inp = branch.astype (float)
branch_inp [:,2] = branch_inp [:,2] + branch_inp [:,3] +branch_inp [:,4]
branch_inp =branch_inp [:,:3] 
branch_inp [branch_inp [:,0]==1, 1:3] = 0.5
branch_inp [:,0] = 1

kernel0 = gpflow.ekernels.RBF (1, variance=0.1)
kernel1 = gpflow.ekernels.RBF (1, variance=0.1)
kernel2 = gpflow.ekernels.RBF (1, variance=0.1)
BP1 = gpflow.params.Parameter ([1.])
cp_kernel1 = CP.tanh_CP ( [bu.zero_kern (), kernel1], locations=BP1)
cp_kernel2 = CP.tanh_CP ( [bu.zero_kern (), kernel2], locations=BP1)
kern_list = gpflow.ParamList([kernel0, cp_kernel1, cp_kernel2])

# ----------equal assignment to main + zero main----------
# for 'model_dict_3b_tb'
branch_num = np.unique (branch.values)
branch = branch.values.reshape ([-1, 1]) == branch_num.reshape ([1, -1])
branch_inp = branch.astype (float)
branch_inp [branch_inp [:,0]==1, 1:3] = 0.5
branch_inp [:,0] = 1

kernel0 = gpflow.ekernels.RBF (1, variance=0.1)
kernel1 = gpflow.ekernels.RBF (1, variance=0.1)
kernel2 = gpflow.ekernels.RBF (1, variance=0.1)
BP1 = gpflow.params.Parameter ([1.])
cp_kernel1 = CP.tanh_CP ( [kernel0, kernel1], locations=BP1)
cp_kernel2 = CP.tanh_CP ( [kernel0, kernel2], locations=BP1)
kern_list = gpflow.ParamList([bu.zero_kern (), cp_kernel1, cp_kernel2])

# ----------1 RBF 2 CP----------
# 'mode_dict_3b_cp'
branch_num = np.unique (branch.values)
branch = branch.values.reshape ([-1, 1]) == branch_num.reshape ([1, -1])
branch_inp = branch.astype (float)
branch_inp [:,0] = 1

kernel0 = gpflow.ekernels.RBF (1, variance=0.1)
kernel1 = gpflow.ekernels.RBF (1, variance=0.1)
kernel2 = gpflow.ekernels.RBF (1, variance=0.1)
BP1 = gpflow.params.Parameter ([1.])
cp_kernel1 = CP.tanh_CP ( [bu.zero_kern(), kernel1], locations=BP1)
cp_kernel2 = CP.tanh_CP ( [bu.zero_kern(), kernel2], locations=BP1)
kern_list = gpflow.ParamList([kernel0, cp_kernel1, cp_kernel2])

# 'model_dict_3b_untied': no tieing the branch points

# ----------successful choice----------
# 1 RBF, 2 CP kernels tied
# duplicate the main branch, use it to fit two kernels
# file: 'model_dict_expand_norm_prior.hdf5'

# 1 RBF, 2 CP kernels tied steepness and branch point
# positive constraint on steepness
# duplicate the main branch, constrain pseudotime
# file: 'module_dict_norm_tie_steep.hdf5'

# 1 RBF, 2 CW kernels tied steepness and branch point
# no constraint on steepness
# duplicate the main branch, constrain pseudotime
# file: 'module_dict_CW_tie_all.hdf5'
