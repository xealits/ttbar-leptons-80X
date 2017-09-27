#from job_handling import *
from job_handling import fetch_n_store_dset, fetch_dset_presence
import os
import argparse
import json, yaml
import logging



parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "Store info on dsets files localy. Default action is to update presence only.",
    epilog = "Example:\n$ python store_dsets_info.py ./dsets/ bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json"
    )

parser.add_argument("dsets_dir",   help="the local directory for dsets files info")
parser.add_argument("dsets",       help="the filename (relational path) or the DAS name of the dsets yaml or json with dtag-dset targets for jobs")
parser.add_argument("--all",       help="update presence of a dset and its' files",
        action = "store_true")
parser.add_argument("--dset",      help="treat input as DAS name",
        action = "store_true")


args = parser.parse_args()


# Prepare proxy file for DAS

if 'X509_USER_PROXY' not in os.environ:
    proxy_file = "./x509_proxy"
    logging.info("Initializing proxy in " + proxy_file)
    #status, output = commands.getstatusoutput('voms-proxy-init --voms cms') # instead of voms-proxy-init --voms cms --noregen
    # os.system handles the hidden input well
    # subprocess.check_output(shell=True) and .Popen(shell=True)
    # had some issues
    status = os.system('voms-proxy-init --voms cms --out ' + proxy_file)

    if status !=0:
        logging.error("Proxy initialization failed:")
        logging.error("    status = " + str(status))
        logging.error("    output =\n" + output)
        exit(2)

    #os.system('export X509_USER_PROXY=' + proxy_file)
    #my_env = os.environ.copy()
    #my_env["X509_USER_PROXY"] = proxy_file
    os.environ["X509_USER_PROXY"] = os.path.abspath(proxy_file)

    logging.info("Proxy is ready.")
else:
    logging.info("using proxy at X509_USER_PROXY")


storing_action = fetch_n_store_dset if args.all else fetch_dset_presence

if args.dset:
    # the dataset DAS name is given on command line
    logging.info("Storing dset %s" % args.dsets)
    storing_action(args.dsets_dir, args.dsets)

elif ".json" in args.dsets:
    # load the dsets from the json file
    with open(args.dsets) as d_f:
        dsets = json.load(d_f)

    logging.info("To local directory %s" % args.dsets_dir)
    for dset in [d['dset'] for dset_group in dsets['proc'] for d in dset_group['data']]:
        logging.info("Storing dset %s" % dset)
        storing_action(args.dsets_dir, dset)

else:
    # load yaml file
    # which has different structure
    with open(args.dsets) as d_f:
        dsets = yaml.load(d_f)

    logging.info("To local directory %s" % args.dsets_dir)
    #for dset in [d['dset'] for dset_group in dsets['proc'] for d in dset_group['data']]:
    for dset in [d for dset_group in dsets for d in dset_group['dsets']]:
        logging.info("Storing dset %s" % dset)
        storing_action(args.dsets_dir, dset)

