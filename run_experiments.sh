
echo "Useful commands:
wc -l experiments/*/err.txt
wc -l experiments/*/out.txt
"
conda activate cell_type_grn_transfer
for experiment in `ls experiments -1` 
do
    echo "Starting ${experiment}"
    # python -m memory_profiler 
    # mprof run --include-children --multiprocess
    python src/do_one_experiment.py --experiment_name $experiment --amount_to_do missing_models --save_trainset_predictions \
        > experiments/$experiment/out.txt 2> experiments/$experiment/err.txt 
done