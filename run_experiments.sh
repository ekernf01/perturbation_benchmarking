source "${HOME}/mambaforge/etc/profile.d/conda.sh"
conda activate ggrn
for experiment in 1.4.2_7 1.4.2_8 1.4.2_9 1.4.2_10 1.4.2_11 1.4.2_12
do
    echo "Starting ${experiment}"
    echo "Monitor progress:
less experiments/${experiment}/err.txt
less experiments/${experiment}/stdout.txt
"
    python do_one_experiment.py --experiment_name $experiment --amount_to_do missing_models --save_trainset_predictions \
        > experiments/$experiment/stdout.txt 2> experiments/$experiment/err.txt

    # To move results back to EK laptop via s3
    echo "aws s3 sync experiments/${experiment} s3://cahanlab/eric.kernfeld/eric_laptop/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/experiments/${experiment}"
    echo "aws s3 sync s3://cahanlab/eric.kernfeld/eric_laptop/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/experiments/${experiment} experiments/${experiment}"
done
