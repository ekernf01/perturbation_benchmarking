# I use these steps to build and publish the Docker image for this project. 
docker build -t ekernf01/perturbation_benchmarking:v1 .
docker run --name perturbation_benchmarking -p 80:80 -d ekernf01/perturbation_benchmarking:v1
docker images
docker tag ekernf01/perturbation_benchmarking:v1 ekernf01/perturbation_benchmarking:v1-release
docker push ekernf01/perturbation_benchmarking:v1-release
