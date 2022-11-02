
echo "Useful commands:
wc -l experiments/*/err.txt
wc -l experiments/*/out.txt
"
conda activate cell_type_grn_transfer
for experiment in `ls experiments -1` 
do
    echo "Starting ${experiment}"
    python -m \
        cProfile -o experiments/$experiment/profiling.txt \
        src/experimenter.py --experiment_name $experiment --amount_to_do models \
        > experiments/$experiment/out.txt 2> experiments/$experiment/err.txt
done
