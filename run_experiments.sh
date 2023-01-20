
echo "Useful commands:
wc -l experiments/*/err.txt
wc -l experiments/*/out.txt
"
conda activate ggrn
for experiment in 1.8.4_0 #`ls experiments -1` 
do
    echo "Starting ${experiment}"
    python src/do_one_experiment.py --experiment_name $experiment --test_mode --amount_to_do missing_models --save_trainset_predictions \
        > experiments/$experiment/out.txt 2> experiments/$experiment/err.txt

done