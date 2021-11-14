# ==========Obtain BRGPLVM prediction==========
# need to run 1GPLVM.R and 2STREAM.R first
maxiter=${1:-1}
mode=${2:-'training'}

root_dir='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root=$root_dir/'results/'
save_dir=$root/'manuscript/figure2/'
script_root=$root_dir/'TBdev/inst/python/BRGPLVM/'
brgp_dir=$save_dir/'BRGP'
mkdir --parent $brgp_dir

#conda activate ptime
cd $script_root
python branching_tb_R.py --root=$save_dir --maxiter=$maxiter --mode=$mode
