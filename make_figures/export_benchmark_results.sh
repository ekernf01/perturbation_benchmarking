
mkdir ../../evaluation_results
for experiment in `ls ../experiments`
do
    mkdir -p ../../evaluation_results/experiments/${experiment}/outputs
    for path_to_copy in metadata.json outputs/conditions.csv outputs/evaluationPerPert.parquet outputs/evaluationPerTarget.parquet outputs/train_resources outputs/train_walltimes
    do
        cp -r ../experiments/${experiment}/${path_to_copy} ../../evaluation_results/experiments/${experiment}/${path_to_copy}
    done
done
cd .. && python gather_experiment_metadata.py
cp ../all_experiments.tsv ../../evaluation_results/all_experiments.tsv