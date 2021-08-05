#!/bin/bash

# Description: Print script usage message
# Arguments: None
# Returns: None
script_usage() {
  echo "
[USAGE]: This script kicks off a script that merges fusion detection results \
from many tools.

-h * [None] Print this help message.
-t * [top_dir1,top_dir2,...] Top directory names under /igm/projects where the \
script will be kicked off.
-s * [samples_file1,samples_file2,...] Files under top_dir1,top_dir2,... that \
contain the lists of samples to kickoff scripts for. The number of samples \
files provided must be equal to the number of top directories, and they must \
be listed in the same respective order. [Default: samples,samples,...].
"
}
