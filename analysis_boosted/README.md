# Way to setup the codes:
- install Madgraph, then install pythia, delphes insides it
- link the `classes` folder of delphes here (use `ln -s XXXX/classes`)
- modify `run.sh` and replace the `<PATH_TO_DELPHES>` to what it should be
- modify `<PATH_TO_DELPHES_OUTPUT>` in `run_detector_analysis_dump.py` to what it should be, and ensure the `seed` is correctly match the naming
Then you are basically ready to run the script locally (`run_detector_analysis_dump.py`)

# submit condor jobs to make it run faster
- create `logs` file here
- check the `seed` setup, e.g. if you have 100 jobs for the generation, the analysis should also be 100 jobs.
- run `condor_submit submit.job`

# analyze the results
- check `run_detector_analysis_postpro.py` and `print_cutflow.py`. They are example on load and process the analysis output files. It is also easy to process using simple python if you do not need to re-generate the QE observables (in the analysis results, the QE observables already calculated. But raw inputs also saved and could be used for advanced study.)
- note: thet are not aim to be directly run. Modify by yourself.