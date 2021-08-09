#!/bin/bash

# Description: Print script usage message
# Globals: (input) TOOL
# Arguments: None
# Returns: None
script_usage() {
  echo "
[USAGE]: This script generates scripts that contain commands to run $TOOL.

-h * [None] Print this help message and exit.
-t * [top_dir1,top_dir2,...] (Required) Top directory names under \$root_dir \
where the scripts will be generated.
-s * [samples_file1,samples_file2,...] (Optional) Files under top_dir1,\
top_dir2,... that contain the lists of samples to generate scripts for. The \
number of samples files provided must be equal to the number of top \
directories, and they must be listed in the same respective order. [Default: \
samples,samples,...].
-r * [root_dir] (Optional) Root directory under which the top directories are \
located. [Default: /igm/projects]
-q * [sge_queue] (Optional) SGE queue to submit jobs. The fusion callers can \
take a long time to run. It can save money to submit long-running jobs to the \
same queue rather than spread the jobs over many nodes. [Default: all.q].
"
}

# Description: Process command line inputs for a single tool gen script
# Globals: (input) TOOL, (output) TOP_DIRS_ARRAY, (output) SAMPLES_FILES_ARRAY,
#   (output) ROOT_DIR, (output) SGE_QUEUE
# Arguments: See script_usage for details
# Returns: None
get_command_line_input() {
  if [ "$#" = "0" ]; then
    script_usage
    exit 1
  fi

  while getopts ":t:s:r:h" opt; do
    case $opt in
      h)
        script_usage
        exit 0
        ;;
      t)
        local top_dirs
        top_dirs="$OPTARG"
        ;;
      s)
        local samples_files
        samples_files="$OPTARG"
        ;;
      r)
        ROOT_DIR="$OPTARG"
        ;;
      :)
        echo "[ERROR]: Flag -$OPTARG requires an argument. See -h for details."
        exit 1
        ;;
      ?)
        echo "[ERROR]: Invalid option -$OPTARG. See -h for valid options."
        exit 1
        ;;
    esac
  done

  if [[ -z "$top_dirs" ]]; then
    echo "[ERROR]: Please specify a top directory. See -h for details."
    exit 1
  fi

  TOP_DIRS_ARRAY=($(echo $top_dirs | sed "s/,/ /g"))
  local num_top_dirs
  num_top_dirs=${#TOP_DIRS_ARRAY[@]}

  if [[ -z "$samples_files" ]]; then
    for (( i=0; i<${num_top_dirs}; i++ )); do
      SAMPLES_FILES_ARRAY[$i]=samples
    done
  else
    SAMPLES_FILES_ARRAY=($(echo $samples_files | sed "s/,/ /g"))
  fi
  local num_samples_files
  num_samples_files=${#SAMPLES_FILES_ARRAY[@]}

  if [[ $num_top_dirs -ne $num_samples_files ]]; then
    echo "[ERROR]: The number of top directories ($num_top_dirs) isn't equal to"
    echo "the number of samples files ($num_samples_files). See -h for details."
    exit 1
  fi

  if [[ -z "$ROOT_DIR" ]]; then
    ROOT_DIR=""
  fi
}
