#!/bin/bash
seed=$1
THISDIR=${PWD}
. /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/setup.sh
export PYTHIA8DATA=<PATH_MADGRAPH>/HEPTools/pythia8/share/Pythia8/xmldoc
export DISPLAY=
mkdir -p run_job_${seed} && cd run_job_${seed}
cp ${THISDIR}/run.cmd .
cp ${THISDIR}/madspin.antitop.card .
cp ${THISDIR}/madspin.top.card .
cp ${THISDIR}/delphes.card .
sed -i "s#<RND>#${seed}#g" run.cmd 
sed -i "s#<RND>#${seed}#g" madspin.antitop.card
sed -i "s#<RND>#${seed}#g" madspin.top.card
<PATH_MADGRAPH>/bin/mg5_aMC run.cmd
cp PROC_loop_sm-no_b_mass_0/Events/run_01/unweighted_events.lhe.gz top.lhe.gz
mv PROC_loop_sm-no_b_mass_0/Events/run_01_decayed_1/tag_1_delphes_events.root top.reco.root
cp PROC_loop_sm-no_b_mass_0/Events/run_01_decayed_1/unweighted_events.root top.parton.root
cp PROC_loop_sm-no_b_mass_1/Events/run_01/unweighted_events.lhe.gz antitop.lhe.gz
mv PROC_loop_sm-no_b_mass_1/Events/run_01_decayed_1/tag_1_delphes_events.root antitop.reco.root
cp PROC_loop_sm-no_b_mass_1/Events/run_01_decayed_1/unweighted_events.root antitop.parton.root
cp PROC_loop_sm-no_b_mass_0/crossx.html top.crossx.html
cp PROC_loop_sm-no_b_mass_1/crossx.html antitop.crossx.html
tar czf outputs.tgz PROC*/ && rm -rf PROC*/
echo "finished ${seed}"
