conda activate ggrn
for experiment in 1.4.2_2 1.4.2_3 1.4.2_4 # 1.2.2_1 1.2.2_2 1.2.2_3 
do
    echo "Starting ${experiment}"
    echo "Monitor progress:
less experiments/${experiment}/err.txt
less experiments/${experiment}/stdout.txt
"
    python do_one_experiment.py --experiment_name $experiment --amount_to_do missing_models --save_trainset_predictions \
        > experiments/$experiment/stdout.txt 2> experiments/$experiment/err.txt
done