#!/bin/bash

#database connection worked!

# Description: Print script usage message
# Arguments: None
# Returns: None
script_usage() {
  echo "
[USAGE]: This script kicks off a script that merges fusion detection results \
from many tools.

-h * [None] Print this help message.
-o * output directory
-s * sample name
-p * patient id
-f * frequency
-u * [false/true] upload to db indicator (default: false)
"
}

# Description: Process command line inputs
# Arguments: See script_usage for details
# Globals: output_location, sample_name, patient_id, upload_to_db
# Returns: None
get_command_line_input() {
  if [ "$#" = "0" ]; then
    script_usage
    exit 1
  fi

  while getopts ":o:s:p:f:u:h" opt; do
    case $opt in
      h)
        script_usage
        exit 0
        ;;
      o)
        output_location="$OPTARG"
        ;;
      s)
        sample_name="$OPTARG"
        ;;
      p)
        patient_id="$OPTARG"
        ;;
      f)
        frequency="$OPTARG"
        ;;
      u)
        upload_to_db="$OPTARG"
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

  if [[ -z "$output_location" ]]; then
    echo "[ERROR]: Please specify an output location. See -h for details."
    exit 1
  fi

  if [[ -z "$sample_name" ]]; then
    echo "[ERROR]: Please specify a sample name. See -h for details."
    exit 1
  fi

  if [[ -z "$patient_id" ]]; then
    echo "[ERROR]: Please specify a patient id. See -h for details."
    exit 1
  fi
  
  if [[ -z "$upload_to_db" ]]; then
    upload_to_db=false
  elif [[ "$upload_to_db" != "true" && "$upload_to_db" != "false" ]]; then
    echo "[ERROR]: -u must be set to 'true' or 'false'. See -h for details."
    exit 1
  fi

}

# Description: Run the overlap script
# Arguments: None
# Globals: script_dir, output_location, sample_name, patient_id, upload_to_db
# Returns: None
run_overlap() {
  echo $output_location
  echo $sample_name
  working_dir=/working_dir
  mkdir -p $working_dir
  cd $working_dir
  echo $working_dir
  cp /$output_location/fusioncatcher/final-list_candidate-fusion-genes.txt .
  cp /$output_location/fusionmap/results/FusionDetection.FusionReport.Table.txt .
  cp /$output_location/jaffa/jaffa_results.csv .
  cp /$output_location/mapsplice/fusions_well_annotated.txt .
  cp /$output_location/starfusion/star-fusion.fusion_predictions.abridged.tsv .
  cp /$output_location/arriba/fusions.tsv .
  cp /$output_location/cicero/annotated.fusion.txt .
  # set ENVIRONMENT for the R session (used to switch between test/prod DB)
  echo "ENVIRONMENT=$ENVIRONMENT" > .Renviron
  if [[ "$upload_to_db" == "true" ]]; then
    $script_dir/R/upload_fusion_results.R --subject $patient_id --sample $sample_name
  fi
  if [[ -z "$frequency" ]]; then
    $script_dir/R/assemble_results.R --sample $sample_name --outReport $working_dir/overlap_$sample_name.tsv
  else 
  $script_dir/R/assemble_results.R --sample $sample_name --outReport $working_dir/overlap_$sample_name.tsv --frequency $frequency
  fi 
#   cp $working_dir/all_fusions_overlap_$sample_name.tsv /$output_location/all_fusions_overlap_$sample_name.tsv
   cp $working_dir/overlap_$sample_name.tsv /$output_location/overlap_$sample_name.tsv
   cp $working_dir/Singleton_KnownFusions_$sample_name.tsv /$output_location/Singleton_KnownFusions_$sample_name.tsv
   cp $working_dir/filtered_overlap_knownfusionlist_3callers_$sample_name.tsv /$output_location/filtered_overlap_knownfusionlist_3callers_$sample_name.tsv
   cp $working_dir/collapse_filtered_overlap_knownfusionlist_3callers_$sample_name.tsv /$output_location/collapse_filtered_overlap_knownfusionlist_3callers_$sample_name.tsv
}

main() {
  script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
  get_command_line_input "$@"
  run_overlap
}

main "$@"
