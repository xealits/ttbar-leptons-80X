"""
The script submitting jobs to LSF on lxplus.
Takes as input:
 * executable
 * cfg template
 * dsets.json

Creates the job-farm directory.
And launches the jobs.

The job is:

    exec cfg.py

The `cfg.py` is compiled from `cfg template` and `dsets.json`.

The job is submitted to LSF as a shell script `job.sh`.
(Which contains `exec cfg.py` inside somewhere.)
The submittion to the LSF job queue is done from shell with `bsub` command:

    $ bsub ... job.sh


# More info

A job.sh looks like:
    #! /bin/sh
    pwd
    export SCRAM_ARCH=slc6_amd64_gcc493
    export BUILD_ARCH=slc6_amd64_gcc493
    export VO_CMS_SW_DIR=/nfs/soft/cms
    cd /afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_5/src/UserCode/llvv_fwk
    eval `scramv1 runtime -sh`
    cd -
    ulimit -c 0;
    ttbarleps80_eventSelection /afs/cern.ch/user/o/otoldaie/work/private/16/CMSSW_8_0_5/src/UserCode/llvv_fwk/data/output/ttbar-leps-80X_1//MC13TeV_TTJets_powheg_514_cfg.py


Job submittion is done with:
    bsub -q 8nh -R "pool>30000" -J ttbarleps80_eventSelection_MC13TeV_TTJets_powheg0000_ -oo /afs/cern.ch/user/o/otoldaie/work/private/16/CMSSW_8_0_5/src/UserCode/llvv_fwk/data/output/ttbar-leps-80X_1//FARM/outputs/0000_ttbarleps80_eventSelection_MC13TeV_TTJets_powheg.cout '/afs/cern.ch/user/o/otoldaie/work/private/16/CMSSW_8_0_5/src/UserCode/llvv_fwk/data/output/ttbar-leps-80X_1//FARM/inputs/job0000_ttbarleps80_eventSelection_MC13TeV_TTJets_powheg.sh'


The cfg.template has the fields `@field` to fill in:
    theLumiMask = path.expandvars(@lumiMask)
    
    runProcess = cms.PSet(
        dtag  = cms.string("@dtag"),
        job_num  = cms.string("@job_num"),
        input = cms.untracked.vstring("@input"),
        outfile = cms.string("@outfile"),

--- TODO: make a pythonic template, with {} instead of the weird @foo.

The dsets.json file structured as:
    {
      "proc":[
    
        {
          "tag":"data",
          "isdata":true,
          "color":1,
          "fill":0,
          "comment":"real split 17",
           "data":[
    
            { "dtag":"Data13TeV_SingleMuon2016B_PromptRecoV2", "split_NoUse_TestingAutomaticOne":50, "xsec":1.0, "br":[ 1.0 ],
              "dset":"/SingleElectron/Run2016B-PromptReco-v2/MINIAOD",
              "lumiMask":"/afs/cern.ch/user/o/otoldaie/work/private/project/CMSSW_7_6_3/src/UserCode/llvv_fwk/data/json/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt" }
    
          ]
        }
    
      ]
    }


"""

#TODO: import argparse

from sys import argv, exit
import json, os
import commands # for some reason subprocess.check_output does not work with das_client
import shutil



if len(argv) != 5:
    print("Usage:\njob_submit.py executable_filename cfg.py_template_filename dsets.json outdirname")
    exit(1)


# Process input

exec_name, cfg_templ_filename, dsets_json_filename, outdirname = argv[1:]
outdirname = os.path.abspath(outdirname.strip()) + '/' # just in case

with open(cfg_templ_filename) as t_f, open(dsets_json_filename) as d_f:
    dsets = json.load(d_f)
    cfg_templ = t_f.read()


job_dir = outdirname + "/FARM/inputs/"
job_outs = outdirname + "/FARM/outputs/"


# Make job FARM directories

print("Making the task (job batch) area in {}".format(outdirname))

try:
    os.makedirs(job_dir)
    os.makedirs(job_outs)
except OSError, e:
    if e.errno == 17:
        print("The directories exist, no need to make them.")
    else:
        raise

project_dir = os.getcwd() + '/' # just in case
print("Project dir  is set to:", project_dir)

# Job templates

# TODO: add SCRAM_ARCH and other parameters
job_template = """#!/bin/sh
pwd
export SCRAM_ARCH={SCRAM_ARCH}
export BUILD_ARCH={SCRAM_ARCH}
export VO_CMS_SW_DIR=/nfs/soft/cms
cd {{project_dir}}
eval `scramv1 runtime -sh`
cd -
ulimit -c 0;
{{exec_name}} {{job_cfg}}
""".format(**os.environ)

job_command_template = """bsub -q 8nh -R "pool>30000" -J {job_name} -oo {job_stdout} '{jobsh}'"""

# could use list comprehension instead of this
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

# Create proxy for working with LSF jobs
# Go through dsets.json
# Create jobs.sh
# Submit them

proxy_file = outdirname + "/FARM/inputs/x509_proxy"
print("Initializing proxy in " + proxy_file)
#status, output = commands.getstatusoutput('voms-proxy-init --voms cms') # instead of voms-proxy-init --voms cms --noregen
# os.system handles the hidden input well
# subprocess.check_output(shell=True) and .Popen(shell=True)
# had some issues
status = os.system('voms-proxy-init --voms cms --out ' + proxy_file)

if status !=0:
    print("Proxy initialization failed:")
    print("    status = " + str(status))
    print("    output =\n" + output)
    exit(2)

os.system('export X509_USER_PROXY=' + proxy_file)

print("Proxy is ready.")


# finding local computing resourses
local_tier = "foo"
hostname = commands.getstatusoutput("hostname -f")[1]
if "cern.ch" in hostname: local_tier = "T2_CH_CERN"

# default file server is the global server:
file_server = "root://cms-xrd-global.cern.ch/"

n_files_per_job = 5

for dset_group in dsets['proc']:
    isdata = dset_group.get('isdata', False)
    # TODO: other parameters?
        # two required parameters of a dset
    for d in dset_group['data']:
        dtag = d['dtag']
        dset = d['dset']
        print("Submitting dtag " + dtag)
        print(dset)

        lumiMask = d.get('lumiMask', '')
        # TODO: other parameters as well?

        # get files of the dset
        status, out = commands.getstatusoutput('das_client --limit=0 --query="file dataset={}"'.format(dset))
        #out_rows = out.split('\n')
        dset_files = out.split('\n')

        if status != 0 or dset_files < 1:
            print("Failed fetch _files_ of {} dataset".format(dset))
            print("das_client status = " + str(status))
            print("           output = " + out)
            print("Continue to other dsets")
            continue

        #dset_files = out_rows[3:]
        print("Found {} files. Splitting {} per job.".format(len(dset_files), n_files_per_job))

        # Finding full local sample
        status, out = commands.getstatusoutput('das_client --query="site dataset={}" --format=JSON '.format(dset))
        #sites = out.split('\n')
        if status != 0 or len(out) < 1:
            print("Failed to fetch _sites_ of {} dataset".format(dset))
            print("das_client status = " + str(status))
            print("           output = " + out)
            print("Continue to other dsets")
            continue
        sites_info = json.loads(out)
        sites = []
        for i in sites_info['data']:
            print(i['site'])
            sites.append(i['site'][0]['name'] + ' ' + i['site'][0]['dataset_fraction'])
        #sites = [i['site'][0]['name'] + ' ' + i['site'][0]['dataset_fraction'] for i in 

        if any(local_tier in s and "100.00%" in s for s in sites):
            print("The full dataset is found on local tier " + local_tier)
            file_server = "root://eoscms//eos/cms/"
            print("using " + file_server + " fileserver")
        else:
            file_server = "root://cms-xrd-global.cern.ch/"
            print("using default " + file_server + " fileserver")

        #sites = sites_rows[3:]

        # the list of bsubs job submition commands
        dset_bsubs = []

        # make jobs per file chunks
        for i, job_chunk in enumerate(chunks(dset_files, n_files_per_job)):
            job_name = exec_name + '_' + dtag + '_' + str(i) + '_'

            # if the site of dataset is CERN -- can run directly on from eos
            #root://eoscms//eos/cms/
            # finding "physical" files via root://cms-xrd-global.cern.ch/
            input_files = ',\n'.join('"%s/%s"' % (file_server, s) for s in job_chunk)

            #if i < 5: print(input_files)

            job_cfg = cfg_templ.format(input = input_files, lumiMask = lumiMask, dtag = dtag, outdir = outdirname, job_num = i, isMC = not isdata, outfile = outdirname + dtag + '_' + str(i) + '.root', project_dir = project_dir)
            # job_cfg = cfg_templ.format(job_files = job_chunk, lumiCert = lumiMask, dtag = dtag, outdir = outdirname, jobID = i, isdata = isdata)

            job_cfg_filename = outdirname + dtag + '_' + str(i) + '_cfg.py'
            with open(job_cfg_filename, 'w') as cfg_file:
                cfg_file.write(job_cfg)

            job_sh = job_template.format(project_dir = project_dir, exec_name = exec_name, job_cfg = job_cfg_filename)
            job_sh_filename = outdirname + '/FARM/inputs/' + dtag + '_' + str(i) + '.sh'
            with open(job_sh_filename, 'w') as sh_file:
                sh_file.write(job_sh)
                os.system("chmod 777 " + sh_file.name)

            job_stdout = outdirname + '/FARM/outputs/' + dtag + '_' + str(i) + '.cout'
            dset_bsubs.append(job_command_template.format(job_name = job_name, job_stdout = job_stdout, jobsh = job_sh_filename))

        print("Created {} jobs. Ready to submit.".format(len(dset_bsubs)))
        with open(outdirname + '/FARM/inputs/' + dtag + '_' + 'bsub.sh', 'w') as bsub_file:
            bsub_file.write('\n'.join(dset_bsubs))
            print("Wrote the bsub commands to {}".format(bsub_file.name))

        print("Ready to submit.")
        #print("No submit in the test run.")
        #print("printing the bsubs instead")

        for bsub in dset_bsubs:
            #print(bsub)
            print(commands.getstatusoutput(bsub))








