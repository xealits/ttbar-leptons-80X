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
        outdir = cms.string("@outdir"),

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

import argparse
import logging

from sys import argv, exit
import json
import os
from datetime import datetime, timedelta

import commands # for some reason subprocess.check_output does not work with das_client
import subprocess # trying it again
import shutil


logging.basicConfig(level=logging.DEBUG)

# Job templates
hostname = commands.getstatusoutput("hostname -f")[1]
logging.debug("hostname = %s" % hostname)
hostname = '.'.join(hostname.split('.')[1:])
logging.debug("hostname = %s" % hostname)

site_cfgs = { 'cern.ch': {'proxy_filename': '/afs/cern.ch/user/o/otoldaie/x509_proxy', 'VO_CMS_SW_DIR': '/nfs/soft/cms',
                          'job_template': """#!/bin/sh
pwd
export X509_USER_PROXY={x509_proxy}
export SCRAM_ARCH={SCRAM_ARCH}
export BUILD_ARCH={SCRAM_ARCH}
export VO_CMS_SW_DIR={VO_CMS_SW_DIR}
cd {{project_dir}}
eval `scramv1 runtime -sh`
cd -
ulimit -c 0;
{{exec_name}} {{job_cfg}}
""",
                          'job_bsub_template': """bsub -q 8nh -R "pool>30000" -J {job_name} -oo {job_stdout} '{jobsh}'"""},
        'ncg.ingrid.pt': {'proxy_filename': '/exper-sw/cmst3/cmssw/users/olek/x509_proxy', 'VO_CMS_SW_DIR': '/cvmfs/cms.cern.ch',
                          'job_template': """#!/bin/sh
pwd
export X509_USER_PROXY={x509_proxy}
export SCRAM_ARCH={SCRAM_ARCH}
export BUILD_ARCH={SCRAM_ARCH}
export VO_CMS_SW_DIR={VO_CMS_SW_DIR}
source $VO_CMS_SW_DIR/cmsset_default.sh
export CMS_PATH=$VO_CMS_SW_DIR
cd {{project_dir}}
eval `scramv1 runtime -sh`
cd -
ulimit -c 0;
{{exec_name}} {{job_cfg}}
""",
                          'job_bsub_template': """qsub -l h_vmem=3G '{jobsh}'"""},
        }

vars_for_the_job = dict(os.environ)
vars_for_the_job.update({'x509_proxy':site_cfgs[hostname]['proxy_filename'], 'VO_CMS_SW_DIR':site_cfgs[hostname]['VO_CMS_SW_DIR']})
job_template = site_cfgs[hostname]['job_template'].format(**vars_for_the_job)

job_command_template = site_cfgs[hostname]['job_bsub_template']


def restore_dset_files_info(dsets_dir, dset):
    """restore_dset_files_info(dsets_dir, dset)

    format:
    dsets_dir/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8,RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1,MINIAODSIM/{files,presence/{date}/[file_servers,block_servers]}

    returns [files], [file_servers]
        --- [file_servers] contains names of tiers (T2_CH_CERN etc) with 100\% of the files (blocks?)

    EXAMPLE:
    restore_dset_files_info("./dsets/", "/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD")
    """
    # t.strftime("%d-%m-%Y")

    dset_dir = dsets_dir + '/' + dset[1:].replace('/', ',')
    logging.info(dset_dir)

    # get files of the dataset
    with open(dset_dir + '/files', 'r') as fs:
        files = [f.strip() for f in fs.readlines()]

    # get the newest directory on presence of files in tiers
    # (tiers with of 100% files -- try 100% of blocks)
    all_subdirs = [dset_dir + '/presence/' + d for d in os.listdir(dset_dir + '/presence') if os.path.isdir(dset_dir + '/presence/' + d) and d[0] != '.']
    logging.debug(all_subdirs)
    newest_dir = max(all_subdirs, key=os.path.getmtime)
    logging.info(newest_dir)

    stored_dset_time = datetime.strptime(os.path.basename(newest_dir), "%d-%m-%Y")
    #datetime.datetime.now() - datetime.timedelta(days=7) < datetime.datetime.strptime("15-11-2016", "%d-%m-%Y")
    if datetime.now() - timedelta(days=7) > stored_dset_time:
        logging.info("restoring a more than a week old dset")
        # TODO: spawn a process to update the record

    #file_server = "root://cms-xrd-global.cern.ch/" # default
    #if os.path.isfile(newest_dir + '/T2_CH_CERN'):
        #file_server = "root://eoscms//eos/cms/"
    file_servers = [d for d in os.listdir(newest_dir + '/file_servers') if os.path.isfile(newest_dir + '/file_servers/' + d) and d[0] != '.']
    logging.debug("100p file_servers: %s" % str(file_servers))

    return files, file_servers

def check_if_stored_dset_info_today(dsets_dir, dset):
    today = datetime.strftime(datetime.now(), "%d-%m-%Y")

    dset_dir = dsets_dir + '/' + dset[1:].replace('/', ',') + '/presence/' + today

    if os.path.isdir(dset_dir):
        logging.debug("Already have stored this dataset info today")
        return True

    return False

def fetch_dset_presence(dsets_dir, dset):
    """fetch_n_store_dset(dsets_dir, dset)
    fetch_n_store_dset("./dsets/", "/SingleElectron/Run2016D-23Sep2016-v1/MINIAOD")
    """
    if check_if_stored_dset_info_today(dsets_dir, dset):
        return True

    # fetch info from DAS
    dset_file_servers = get_dset_presence(dset)

    #if not (dset_files and dset_file_servers):
        #print("Failed to fetch dset files info from DAS")
        #return False

    return store_dset_presence(dsets_dir, dset, dset_file_servers)


def fetch_n_store_dset(dsets_dir, dset):
    """fetch_n_store_dset(dsets_dir, dset)
    fetch_n_store_dset("./dsets/", "/SingleElectron/Run2016D-23Sep2016-v1/MINIAOD")
    """
    if check_if_stored_dset_info_today(dsets_dir, dset):
        return True

    # fetch info from DAS
    dset_files, dset_file_servers = get_dset_files(dset), get_dset_presence(dset)

    return store_dset_files(dsets_dir, dset, dset_files) and store_dset_presence(dsets_dir, dset, dset_file_servers)


def store_dset_files(dsets_dir, dset, dset_files):
    """store_dset_presence(dsets_dir, dset, dset_files)

    format:
    dsets_dir/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8,RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1,MINIAODSIM/files

    returns True # on success

    EXAMPLE:
    store_dset_files_info("./dsets/", "/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD", ["file1", "file2"])
    """

    dset_dir = dsets_dir + '/' + dset[1:].replace('/', ',')
    if not os.path.exists(dset_dir):
        os.makedirs(dset_dir)

    with open(dset_dir + '/files', 'w') as f:
            f.write('\n'.join(dset_files))

    return True


def store_dset_presence(dsets_dir, dset, dset_file_servers):
    """store_dset_presence(dsets_dir, dset, dset_file_servers)

    format:
    dsets_dir/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8,RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1,MINIAODSIM/presence/15-11-2016/file_servers/T2_CH_CERN

    returns True # on success

    EXAMPLE:
    store_dset_files_info("./dsets/", "/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD", ["T2_CH_CERN", "T2_DE_DESY"])
    """

    today = datetime.strftime(datetime.now(), "%d-%m-%Y")

    dset_dir = dsets_dir + '/' + dset[1:].replace('/', ',') + '/presence/' + today

    # let's leave it fail here, the check is done by hand or in a script above
    #if os.path.isdir(dset_dir):
        #print("Already have stored this dataset info today")
        #return True
    if not os.path.exists(dset_dir + '/file_servers/'):
        os.makedirs(dset_dir + '/file_servers/')

    # if any 100% full servers are found -- touch the file with server name
    for s in dset_file_servers:
        s_file = open(dset_dir + '/file_servers/' + s, 'a')
        s_file.close()

    return True


def store_dset_files_info(dsets_dir, dset, dset_files, dset_file_servers):
    """store_dset_files_info(dsets_dir, dset, dset_files, dset_file_servers)

    format:
    dsets_dir/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8,RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1,MINIAODSIM/15-11-2016/T2_CH_CERN

    returns True # on success

    EXAMPLE:
    store_dset_files_info("./dsets/", "/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD", ["file1", "file2"], ["T2_CH_CERN", "T2_DE_DESY"])
    """

    today = datetime.strftime(datetime.now(), "%d-%m-%Y")

    dset_dir = dsets_dir + '/' + dset[1:].replace('/', ',') + '/presence/' + today

    # let's leave it fail here, the check is done by hand or in a script above
    #if os.path.isdir(dset_dir):
        #print("Already have stored this dataset info today")
        #return True
    if not os.path.exists(dset_dir + '/file_servers/'):
        os.makedirs(dset_dir + '/file_servers/')

    # if any 100% full servers are found -- touch the file with server name
    for s in dset_file_servers:
        s_file = open(dset_dir + '/file_servers/' + s, 'a')
        s_file.close()

    if not dset_files:
        print("Failed to fetch files.")
        failed_files = open(dset_dir + '/files_FAILED', 'w')
        failed_files.close()
        return False

    with open(dset_dir + '/files', 'w') as f:
            f.write('\n'.join(dset_files))

    return True


def get_dset_presence(dset):
    """get_dset_presence(dset)

    EXAMPLE:
    get_dset_presence("/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD")
    """
    # Get fileserver info (100% of files)

    # Finding full local sample
    logging.debug("Command:")
    logging.debug('./das_client --query="site dataset={}" --format=JSON '.format(dset))
    status, out = commands.getstatusoutput('./das_client --query="site dataset={}" --format=JSON '.format(dset))
    #sites = out.split('\n')
    if status != 0 or len(out) < 1:
        logging.error("Failed to fetch _sites_ of {} dataset".format(dset))
        logging.error("das_client status = " + str(status))
        logging.error("           output = " + out)
        logging.error("Continue to other dsets")
        return None, None
    sites_info = json.loads(out)
    sites = []

    for i in sites_info['data']:
        logging.debug(i['site'])
        #x = i['site'][0]
        #x.update(i['site'][1])
        x = {}
        for s in i['site']:
            x.update(s)
        sites.append((x['name'], x.get('dataset_fraction')))

    logging.info("Found sites and completeness:")
    for s in sites:
        logging.info(s)

    '''
    if any(local_tier in s and "100.00%" in s for s in sites):
        print("The full dataset is found on local tier " + local_tier)
        file_server = "root://eoscms//eos/cms/"
        print("using " + file_server + " fileserver")
    else:
        file_server = "root://cms-xrd-global.cern.ch/"
        print("using default " + file_server + " fileserver")
    '''
    # return all tier names, which have 100% of the dataset
    return [s[0] for s in sites if "100.00%" in s]



def get_dset_files(dset):
    """get_dset_files(dset)

    EXAMPLE:
    get_dset_files("/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD")
    """
    #status, out = commands.getstatusoutput('das_client --limit=0 --query="file dataset={}"'.format(dset))
    # trying os.system to get the environment variable of X509_USER_PROXY propagate to das call:
    #status, out = os.system('das_client --limit=0 --query="file dataset={}"'.format(dset)), '<output is at stdout>'
    # -- no! the output of the command is the list of files
    #status, out = 0, subprocess.check_output( 'export X509_USER_PROXY=' + proxy_file + ' && das_client --limit=0 --query="file dataset={}"'.format(dset), shell=True)
    logging.debug("Command:")
    logging.debug('./das_client --limit=0 --query="file dataset={}"'.format(dset))
    status, out = 0, subprocess.check_output( './das_client --limit=0 --query="file dataset={}"'.format(dset), shell=True)

    # the correct approach
    #p = subprocess.Popen(['das_client', '--limit=0', '--query="file dataset={}"'.format(dset)], env=my_env)
    #out, err = p.communicate()
    #status = 0 # TODO: extract exit status

    # TODO: add the status of the commend exit, check_output raises an Error if the status code is bad
    out_rows = out.split('\n')
    dset_files = out.strip().split('\n')
    
    if status != 0 or dset_files < 1:
        logging.error("Failed fetch _files_ of {} dataset".format(dset))
        logging.error("das_client status = " + str(status))
        logging.error("           output = " + out)
        logging.error("Continue to other dsets")
        return None, None

    return dset_files




#def get_dset_files(dset):
def get_dset(dset):
    """get_dset(dset)

    EXAMPLE:
    get_dset("/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD")
    """
    #status, out = commands.getstatusoutput('das_client --limit=0 --query="file dataset={}"'.format(dset))
    # trying os.system to get the environment variable of X509_USER_PROXY propagate to das call:
    #status, out = os.system('das_client --limit=0 --query="file dataset={}"'.format(dset)), '<output is at stdout>'
    # -- no! the output of the command is the list of files
    #status, out = 0, subprocess.check_output( 'export X509_USER_PROXY=' + proxy_file + ' && das_client --limit=0 --query="file dataset={}"'.format(dset), shell=True)
    print("Command:")
    print('./das_client --limit=0 --query="file dataset={}"'.format(dset))
    status, out = 0, subprocess.check_output( './das_client --limit=0 --query="file dataset={}"'.format(dset), shell=True)

    # the correct approach
    #p = subprocess.Popen(['das_client', '--limit=0', '--query="file dataset={}"'.format(dset)], env=my_env)
    #out, err = p.communicate()
    #status = 0 # TODO: extract exit status

    # TODO: add the status of the commend exit, check_output raises an Error if the status code is bad
    out_rows = out.split('\n')
    dset_files = out.strip().split('\n')
    
    if status != 0 or dset_files < 1:
        print("Failed fetch _files_ of {} dataset".format(dset))
        print("das_client status = " + str(status))
        print("           output = " + out)
        print("Continue to other dsets")
        return None, None
    #dset_files = None # FIXME: remove

    # And now find the fileserver

    #def get_dset_site(dset):
    # Finding full local sample
    print("Command:")
    print('./das_client --query="site dataset={}" --format=JSON '.format(dset))
    status, out = commands.getstatusoutput('./das_client --query="site dataset={}" --format=JSON '.format(dset))
    #sites = out.split('\n')
    if status != 0 or len(out) < 1:
        print("Failed to fetch _sites_ of {} dataset".format(dset))
        print("das_client status = " + str(status))
        print("           output = " + out)
        print("Continue to other dsets")
        return None, None
    sites_info = json.loads(out)
    sites = []
    #print(sites_info['data'])
    for i in sites_info['data']:
        print(i['site'])
        #x = i['site'][0]
        #x.update(i['site'][1])
        x = {}
        for s in i['site']:
            x.update(s)
        sites.append((x['name'], x.get('dataset_fraction')))
    #sites = [i['site'][0]['name'] + ' ' + i['site'][0]['dataset_fraction'] for i in 

    print("Found sites and completeness:")
    for s in sites:
        print(s)

    '''
    if any(local_tier in s and "100.00%" in s for s in sites):
        print("The full dataset is found on local tier " + local_tier)
        file_server = "root://eoscms//eos/cms/"
        print("using " + file_server + " fileserver")
    else:
        file_server = "root://cms-xrd-global.cern.ch/"
        print("using default " + file_server + " fileserver")
    '''
    # return all tier names, which have 100% of the dataset
    return dset_files, [s[0] for s in sites if "100.00%" in s]

known_file_tiers = {"cern.ch": ("T2_CH_CERN", "root://eoscms//eos/cms/")}
"""dictionary holding
    {hostname: (TIER_NAME, tier_file_server), ...}

EXAMPLE:
    known_file_tiers = {"cern.ch": ("T2_CH_CERN", "root://eoscms//eos/cms/")}

Use with default file_server:
    local_tier, local_file_server = known_file_tiers.get(hostname, ("CMS_DEFAULT_GLOBAL_XRD", "root://cms-xrd-global.cern.ch/"))
"""

# default file server is the global server:
default_file_server = "root://cms-xrd-global.cern.ch/"

