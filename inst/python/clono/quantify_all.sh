#! /bin/bash
trial=${1:-'clono'}
# root: the root directory containing both the python venv and data
root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
data_dir=${2:-$root'results/imaging/clono'}
env_dir=${3:-$root'amnioid/cellprofiler'}

for file in $(ls $data_dir/$trial)
do
        echo "quantifying "$file
        ./quantify.sh $file $trial $data_dir $env_dir
done

./combine_all.R $data_dir $trial

: <<'END_OF_DOCS'
Argument:
1: trial, the folder containing the data to quantify in $data_dir
2: data_dir, where all the image data is located. Each image must be stored
separately as .tif file in 2D, that is, just height and width. Every image does
not have to be of the same size.
3: env_dir, path to the python virtual environment containing cellprolifer.

Description: 
Run this script in its current directory
example folder structure:
root
----amnioid/cellprolifer #virtual environment ($env_dir)
----results/imaging #all imaging data ($data_dir)
--------------------clono # trial to quantify ($trial)
-------------------------anxa4 #each condition in the trial
-------------------------tbx3
-------------------------cebpa
------------------------------cebpa_001.tif #tif file containing 2D images
------------------------------cebpa_002.tif 

Output:
A directory will be created as: $data_dir/$trial/proc_data
In this directory, for every condition in the trial, there will be a
'*_sum_.csv' file for the quantification results. This includes the morphology
of every nucleus segmented in this image

Example:
$ root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
$ data_dir=$root'results/imaging'
$ env_dir=$root'amnioid/cellprofiler'
$ ./quantify_all.sh 'clono' $data_dir $env_dir

Output:
END_OF_DOCS
