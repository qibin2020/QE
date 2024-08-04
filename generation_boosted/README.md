# Way to setup the codes:
- install Madgraph, then install pythia, delphes insides it
- modify `setup.sh` and `run.sh`, replace the `<PATH_MADGRAPH>` to what it should be
Then you are basically to run locally (`run.sh`)

# submit condor jobs to make it run faster
- create `logs` file here
- check the `seed` setup, e.g. if you want to run 100 jobs, modify `queue` numbers. Every job is 10k events. This is boosted only.
- run `condor_submit submit.job`
