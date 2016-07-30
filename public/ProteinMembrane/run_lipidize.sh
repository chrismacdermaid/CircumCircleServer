#!/bin/bash
#set -u

## Run lipidize
echo "Running Lipidize"

# added by Anaconda2 4.1.1 installer
export PATH="/home/vagrant/pkg/anaconda2/bin:$PATH" && source activate lipidize

# Check that file exists, sanatize, and execute
[[ -e "${1}" ]] && FNAME=`ls "${1}"` && python lipidize.py ${FNAME} && exit 0

echo "Lipidizei Complete"

exit 1
