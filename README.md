# Establish a trophoblast developmental trajectory

This respository contains the code to reproduce the results in the paper:

![](vignettes/TB_trajectory.png)

## Installation
The package can be installed in R using:
```{r}
devtools::install_github ('Yutong441/TBdev')
```

In case you would like to reuse some of our analysis pipeline, please refer to
the 'vignettes' folder for illustrations.
You may also find the data to reproduce the vignettes and the output [here](https://drive.google.com/drive/folders/1Jz2s33SLmvtXisVPTwNZtDBU4uvEInax?usp=sharing).
To alter the default parameters of the graphics, you may download a copy of the
'data-raw/config_template.R' file, alter it, and load this file into R.

## Process expression matrix
First, the raw expression matrices after alignment are processed using the
scripts in the 'manuscript/clean_data' folder. Each script processes the
expression matrix from one scRNA-seq dataset. An exception is that
'manuscript/clean_data/Stirparo_2018.R' processes three studies: Blakeley_2015,
Petropoulos_2016 and Yan_2013. 

Due to the storage constraint of github, it was not possible to upload all of
our raw data. Our raw data would have been stored in the 'data/' directory. It
contains one folder for each study. Each folder has the raw expression matrix.
The scripts in the 'manuscript/clean_data' folder creates a Seurat object for
each study stored in their respective subfolders. After creating those Seurat
objects, the 'merge_data.R' script integrates all the data and store the final
results in a subfolder in the 'result/' folder.

## Re-annotate the lineages
Next the lineages are re-annotated using the scripts in 'assign_label'. First,
the 'assign_label/assign_preimplant.R' file standardises some of the authors'
labels. For example, the some authors use 'TE' to denote pre-implantation
trophoblasts, while Zhou (2019) uses it to denote trophoblast in general. Those
names are standardised. The overall scheme is detailed in
'assign_label/clean_label.md'. Secondly, the hTSC and hESC data are merged into
the dataset. This is done in the 'add_TSC.R' file. Having cleaned the labels,
we performed clustering analysis in the 'assign_label/change_labels.R' file. 24
clusters are identified and annotated individually and stored in the
'utils/labels/cluster_assignment.csv' file. They are integrated into the
metadata of the merged data.  Because we perform pseudotime analysis in python,
which may not process Seurat object, we store the data as csv. The results from
the pseudotime analysis are loaded into the Seurat object for downstream
analysis.

## Generate the figures in this paper
The scripts used to generate the figures in this paper (except for figure 4)
are in the 'manuscript' folder. Each file creates one figure.

## Pseudotime analysis
The pseudotime analysis workflow requires multiple steps. First, in the root
directory, run the 'assign_label/change_labels.R' script section 'save data for
pseudotime analysis' to obtain the input matrix to GPLVM. It generates 2 files
in the 'data/' folder. One is 'merged_all.csv'. The other is
'merged_meta_all.csv'. Next install tensorflow version 1.13.1 and gpflow
version 1.0.0. In the root directory, run 'inst/python/GPLVM/tb_gplvm.py'. This
script uses GPLVM to infer the 3D latent space for the dataset. Having run this
file, there would be 4 files in the 'result/' folder.  The
'model_no_prior.hdf5' file stores the GPLVM model. The 'PT_no_prior.csv' file
contains the latent space locations for each cell. The 'PT_pred_mean.csv' and
'PT_pred_var.csv' stores the inferred mean and variance for each latent space
location. These data are essential for Elpigraph as implemented in the STREAM
package. Alternatively, you may wish to follow the instruction from
'vignettes/GPLVM.Rmd' to run it in R. However, in the original paper, it was
done in python.

Next, install the python package [STREAM](https://github.com/pinellolab/STREAM). 
After running the file 'inst/python/STREAM/stream_GPLVM.py', 3 files in the
'data/' folder will be generated. The 'stream_branch_labels.csv' contains
branch information. 'stream_pseudotime' contains the pseudotime estimated by
STREAM. Lastly, the 'STREAM_data.csv' file contains scaled expression matrix.
These files are inputs to the B-RGPLVM algorithm.

Lastly, run B-RGPLVM, using the script 'inst/python/BRGPLVM/branching_tb.py'.
Three files will be generated in the 'result/' folder.
'model_dict_matern.hdf5' is the B-RGPLVM model. The 'prediction_matern_500.csv'
contains inferred gene expression from B-RGPLVM. Lastly, the
'infer_pt_matern.csv' contains the pseudotime locations estimated by B-RGPLVM.

However, it is difficult to assemble a pipeline connecting all these three
components. First of all, it is necessary to manually inspect the trajectory
inference results from STREAM. Secondly. GPLVM and B-RGPLVM were implemented in
tensorflow 1. However, the STREAM package depends on tensorflow 2. We had to
create two virtual environments as followed:

```{bash}
# create environment for GPLVM and B-RGPLVM
conda create -n ptime python=3.6
conda activate ptime
conda install -c conda-forge tensorflow==1.13.1
conda install pip
pip install gpflow==1.0.0
conda install -c bioconda anndata
conda deactivate ptime

# create another environment for STREAM
conda create -n scrna python=3.6
conda install -c bioconda STREAM
conda deactivate scrna
```
