#!/bin/bash


# Example to run this pipeline:
# export conda_base="/sc/arion/projects/sealfs01a/stas/conda/local_envs/snakemake3_mamba4"
# export data_path="/sc/arion/projects/sealfs01a/stas/local/snakemake_alternative_splicing/data_to_process"
# export pipeline_emails_to="stas.rirak@mssm.edu, stasrirak.ms@gmail.com"
# # export dry_run=True
# # export help=True
# # export debug=True
# # python2_path="/sc/arion/projects/sealfs01a/stas/conda/local_envs/snakemake3_py2/bin/python"

# ./run_pipeline.sh -h  # to show script help info
# ./run_pipeline.sh     # to run pipeline


# Function to display help text
show_help() {
    echo "Options:"
    echo "  --help, -h    Show this help text"
    echo "No other command-line arguments are accepted. However the following options can be provided through environment variables. All options are optional unless marked as 'required'"
    echo "  conda_base      : (required) path to the base conda environment containing the required version of the snakemake installation for this pipeline"
    echo "  data_path       : (required) path to the folder containing data that will be processed"
    echo "  conda_envs_path : path to the directory containing all additional conda environments (in addition to conda_base) required to run some pipeline steps, i.e. python2 environment. If not provided, this will be set to the parent directory of the conda_base location."
    echo "  snakemake_file  : path to the snakemake file located in the root of the pipeline application. Default value: 'alter_splicing.smk' file located in the same folder as the current script."
    echo "  cluster_config  : path to the snakemake cluster config file to be used with the run. Default value: 'config/cluster_cfg.yaml'. If not provided, the pipeline will create it using default config setting."
    echo "  dry_run         : True/False value. Defines 'dry run' setting of the snakemake pipline. Default value: False."
    echo "  cluster_run     : True/False value. Defines if cluster will be used for running the snakemake pipline. Default value: True."
    echo "  jobs            : Integer value. The number of jobs submitted to the cluster at the same time. Default value: 400."
    echo "  cores_local     : Integer value. Defines number of cores to be used by the pipeline when runs locally (won't affect cluster's settings). Default value set by this script: 40."
    echo "  latency_wait    : Integer value. Defines number of seconds the pipeline will wait for results to be returned from the cluster. Default value set by this script: 15."
    echo "  rerun_incomplete: True/False value. Defines 'rerun-incomplete' argument of the snakemake pipeline. If set, the failed steps from previous runs will be forced to be re-done. Default value: True."
    echo "  keep_going      : True/False value. Defines 'keep-going' argument of the snakemake pipline. If set, it will go on with independent jobs if a job fails. Default value: True."
    echo "  verbose         : True/False value. Defines 'verbose' argument of the snakemake pipline. Default value: True."
    echo "  ml_purge        : True/False value. Defines if the 'ml purge' command will be executed before running the pipeline. Default value: True."
    echo "  help            : True/False value. If set True, snakemake execution will halt and additional help info will be provided by the snakemake file. Default value: False."
    echo "  debug           : True/False value. If set True, 1) bash debug mode will be activated, 2) Snakemake pipeline will preserve all 'temp' objects and would not delete them. Default value: False."
    echo "  additional_args : String value. It can contain any additional snakemake argruments that should be supplied for the snakemake run. Default value is blank."
}

# Check for command-line options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help | -h)
            show_help
            exit 0
            ;;
        *)
            echo "Invalid option: $1"
            show_help
            exit 1
            ;;
    esac
    shift
done

script_dir=$(dirname "$0")  # parent folder path of the current script
project_name=$(basename "$script_dir") # parent folder name 

# snakemake_file="/sc/arion/projects/sealfs01a/stas/local/snakemake_alternative_splicing/alter_splicing.smk"

if [ -z "$snakemake_file" ]; then
    snakemake_file=$script_dir/rrbs.smk
    if [ ! -f "$snakemake_file" ]; then
        echo "Error: The following implicitly identified 'snakemake_file' is not a valid file: $snakemake_file"
    fi
else
    if [ ! -f "$snakemake_file" ]; then
        echo "Error: The following explicitly identified 'snakemake_file' is not a valid file: $snakemake_file"
    fi
fi

snakemake_dir=$(dirname "$snakemake_file")
# orig_dir=$(pwd)
cd $snakemake_dir

# Check if $conda_base is empty or not a valid directory
if [ -z "$conda_base" ] || [ ! -d "$conda_base" ]; then
    echo "Error: 'conda_base' variable is empty or not a valid directory or does not exist. This variable must be set with a path to the conda environment location containing the required version of the snakemake tool to run the current pipeline."
    exit 1
fi

# Check if $data_path is empty or not a valid directory
if [ -z "$data_path" ] || [ ! -d "$data_path" ]; then
    echo "Error: 'data_path' variable is empty or not a valid directory or does not exist. This variable must be set with a path to the location of the project that needs to be processed by the pipeline."
    exit 1
fi

# Check if $conda_envs_path is empty or not a valid directory
if [ -z "$conda_envs_path" ] || [ ! -d "$conda_envs_path" ]; then
    # if not a valid value provided, assign it with the parent directory of the conda_base
    # conda_envs_path=$(dirname "$conda_base")
    conda_envs_dir=$(dirname "$conda_base") # parent folder of the conda_base location
    current_user=$(whoami)
    echo "process runs under user = "$current_user
    conda_envs_path="$conda_envs_dir"/temp_envs/"$project_name"_"$current_user"
    echo "Directory for temp conda environments (conda_envs_path) = "$conda_envs_path
    
    mkdir -p "$conda_envs_path"
    chmod -R g+rwx "$conda_envs_path"
    export conda_envs_path=$conda_envs_path
fi

if [ -z "$cluster_config" ]; then
    # If it's blank, assign a default value
    cluster_config="config/cluster_cfg.yaml"
fi

# Check if $ml_purge is set to false, otherwise execute the 'ml purge' command
if [[ "${ml_purge,,}" == "false" ]]; then  # will convert ml_purge value to lower case 
  echo "skipping ml purge command since the environment variable ml_purge was set to False"
else
  ml purge
  echo "ml purge command was applied"
fi

# Check if $dry_run is set to true
if [[ "${dry_run,,}" == "true" ]]; then  # will convert dry_run value to lower case 
    dry_run_script='-np'
else
    dry_run_script='-p'
fi

# Check if $verbose is false
if [[ "${verbose,,}" == "false" ]]; then  # will convert verbose value to lower case 
    verbose=''
else
    verbose="--verbose"
fi

# check $keep_going
if [[ "${keep_going,,}" == "false" ]]; then  # will convert verbose value to lower case 
    keep_going=''
else
    keep_going="--keep-going"
fi

# If it's blank or not an integer, assign a default value
if [ -z "$jobs" ] || ! [[ "$jobs" =~ ^[0-9]+$ ]]; then
    jobs_param="--jobs 400"
else
    jobs_param="--jobs $jobs"
fi

# If it's blank or not an integer, assign a default value
if [ -z "$cores_local" ] || ! [[ "$cores_local" =~ ^[0-9]+$ ]]; then
    cores_local_param="-c40"
else
    cores_local_param="-c$cores_local"
fi

# If it's blank or not an integer, assign a default value
if [ -z "$latency_wait" ] || ! [[ "$latency_wait" =~ ^[0-9]+$ ]]; then
    latency_wait_param="--latency-wait 15"
else
    latency_wait_param="--latency-wait $latency_wait"
fi

# Check if $rerun_incomplete is set to true
if [ -z "$rerun_incomplete" ] || [[ "${rerun_incomplete,,}" == "true" ]]; then  # will convert the rerun_incomplete value to lower case 
    rerun_incomplete_script='--rerun-incomplete'
else
    rerun_incomplete_script=''
fi

# Check if $additional_args is blank
if [ -z "$additional_args" ]; then  
    additional_snakemake_args=''
else
    additional_snakemake_args="$additional_args"
fi

echo "Starting conda_base environment:  $conda_base"
# start runing pipeline
eval "$(conda shell.bash hook)"  # this line is required to run conda activate (source: https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script)
conda activate $conda_base
echo "conda_base environment was started"

# Check if $debug is false
if [[ "${debug,,}" == "true" ]]; then  # will convert debug value to lower case 
    set -x  # start bash debuging mode
fi

# Check if $cluster_run is false
if [[ "${cluster_run,,}" == "false" ]]; then  # will convert cluster_run value to lower case 
    snakemake $dry_run_script $cores_local_param $latency_wait_param $rerun_incomplete_script --use-conda --conda-prefix $conda_envs_path --directory $data_path -s $snakemake_file $verbose $keep_going $additional_snakemake_args
else
    # original way of running snakemake (without a profile)
    # snakemake $dry_run_script $cores_local_param $latency_wait_param $rerun_incomplete_script --use-conda --conda-prefix $conda_envs_path --directory $data_path -s $snakemake_file $verbose $keep_going $additional_snakemake_args $jobs_param --cluster-config $cluster_config --cluster "bsub -W {cluster.time} -P {cluster.account} -q {cluster.queue} -n {cluster.nCPUs} -R {cluster.resources} -J {cluster.name} -oo {cluster.output} -eo {cluster.error}"
    
    # running snakemake using a local profile
    snakemake $dry_run_script $cores_local_param $latency_wait_param $rerun_incomplete_script --conda-prefix $conda_envs_path --directory $data_path -s $snakemake_file $verbose $keep_going $additional_snakemake_args $jobs_param
fi

if [[ "${debug,,}" == "true" ]]; then  # will convert debug value to lower case 
    set +x  # turn off bash debuging mode
fi

conda deactivate

