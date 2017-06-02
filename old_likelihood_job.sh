#!/bin/sh
pwd
export X509_USER_PROXY=/afs/cern.ch/user/o/otoldaie/x509_proxy
export SCRAM_ARCH=slc6_amd64_gcc530
export BUILD_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/nfs/soft/cms
cd /afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_25/src/UserCode/ttbar-leptons-80X/
eval `scramv1 runtime -sh`
ulimit -c 0;
source old_run_regions.C && source old_likelihood_distrs
