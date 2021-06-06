#! /bin/bash
# run this script in its current directory
root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
condition=${1:-'anxa4'}
trial=${2:-'clono'}

# where the .lif input file is stored and the outputs are directed to
data_dir=${3:-$root'results/imaging/clono'}
# where the python environment containing cellprofier is
env_dir=${4:-$root'amnioid/cellprofiler'}

source $env_dir'/bin/activate'

python3 segment_clono.py $data_dir/$trial/$condition
echo 'finishing segmentation'

python3 get_clone.py $data_dir/$trial/$condition
echo 'finishing nuclei quantification'

python3 summarise.py $data_dir/$trial/ $condition
echo 'finishing combining the data'

./combine_raw.R $data_dir/$trial/ $condition
echo 'finishing cleaning the data'
deactivate
