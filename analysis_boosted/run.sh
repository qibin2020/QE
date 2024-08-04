#!/bin/bash
seed=$1
THISDIR=${PWD}
. /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc12-opt/setup.sh
export DISPLAY=
mkdir -p run_job_${seed} && cd run_job_${seed}
ln -s ../run_detector_analysis_dump.py .
ln -s ../kernel* .
ln -s ../classes/ .
ln -s <PATH_TO_DELPHES>/libDelphes.so
time python run_detector_analysis_dump.py ${seed}
echo "finished ${seed}"
