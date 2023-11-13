source "${HOME}/mambaforge/etc/profile.d/conda.sh"
conda activate ggrn
for experiment in `ls -1 experiments | grep -E "1.4.3_*|1.3.3_*|1.0_*|1.2.2_*|1.4.2_*|1.6.1_*|5_0"`

do
    echo "Starting ${experiment}"
    echo "Monitor progress:
less experiments/${experiment}/err.txt
less experiments/${experiment}/stdout.txt
"
    python do_one_experiment.py --experiment_name $experiment --amount_to_do missing_models \
        > experiments/$experiment/stdout.txt 2> experiments/$experiment/err.txt
done
