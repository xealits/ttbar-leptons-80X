#from job_handling import *
from job_handling import fetch_n_store_dset
import os
import argparse
import json



parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "Store info on dsets files localy.",
    epilog = "Example:\n$ python store_dsets_info.py ./dsets/ bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json"
    )

parser.add_argument("dsets_dir",  help="the local directory for dsets files info")
parser.add_argument("dsets",    help="the filename (relational path) of the dsets json with dtag-dset targets for jobs")

args = parser.parse_args()


# Prepare proxy file for DAS

proxy_file = "./x509_proxy"
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

#os.system('export X509_USER_PROXY=' + proxy_file)
#my_env = os.environ.copy()
#my_env["X509_USER_PROXY"] = proxy_file
os.environ["X509_USER_PROXY"] = os.path.abspath(proxy_file)

print("Proxy is ready.")

'''
for dset_group in args.dsets['proc']:
    for d in dset_group['data']:
        dset = d['dset']
'''

with open(args.dsets) as d_f:
    dsets = json.load(d_f)

print("To local directory %s" % args.dsets_dir)
for dset in [d['dset'] for dset_group in dsets['proc'] for d in dset_group['data']]:
    print("Storing dset %s" % dset)
    fetch_n_store_dset(args.dsets_dir, dset)

