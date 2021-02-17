# run this script in its current directory
root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
# where the .lif input file is stored and the outputs are directed to
data_dir=$root'results/imaging'
# where the python environment containing cellprofier is
env_dir=$root'amnioid/cellprofiler'
condition=${2:-'all'}
channels=${3:-DAPI HLAG CGB TFAP2C}

source $env_dir'/bin/activate'
python3 'lif_to_tiff.py' $data_dir/"$1" "$condition" --c $channels
echo 'finishing converting LIF to TIFF'
python3 segment_cyto.py $data_dir"/data/"$condition
echo 'finishing segmentation'
python3 fluorescence.py $data_dir"/data/"$condition --n $channels
echo 'finishing fluorescence quantification'
python3 summarise.py $data_dir"/data/" $condition
echo 'finishing combining the data'
./combine_raw.R $data_dir $condition
echo 'finishing cleaning the data'
deactivate
