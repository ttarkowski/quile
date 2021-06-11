#!/bin/bash

# Arguments
# 1 - pw.x input file

error_msg="Calculations failed."

# DFT calculations
pw.x < $1 >& ${1}.out

# Cleanup
rm -rf results_${1}/

# Error message and exit in case of some calculations problems
stop=`cat ${1}.out | grep "STOP" | wc -l`
if [ ${stop} -gt 0 ]; then
  echo ${error_msg}
  exit
fi

# The final result: energy or error message
job_done=`cat ${1}.out | grep "JOB DONE" | wc -l`
if [[ ${job_done} == 1 ]]; then
  cat ${1}.out | grep "!" | grep "total energy" | tail -n 1 | awk '{print $5}'
else
  echo ${error_msg}
fi
exit

