source "${HOME}/mambaforge/etc/profile.d/conda.sh"
conda activate ggrn
for experiment in  "1.4.2_1" "1.4.2_2" "1.4.2_3" "1.4.2_5" "1.4.2_6" "1.4.2_7"
#`ls -1 experiments | grep 1.4.2_*` 
do
    echo "Starting ${experiment}"
    echo "Monitor progress:
less experiments/${experiment}/err.txt
less experiments/${experiment}/stdout.txt
"
    python do_one_experiment.py --experiment_name $experiment --amount_to_do missing_models --save_trainset_predictions \
        > experiments/$experiment/stdout.txt 2> experiments/$experiment/err.txt
done
