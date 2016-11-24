# very bad practice:
from job_handling import *


# why __main__ if the whole file is script?
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = "Submit jobs to LSF cluster on lxplus.",
        epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
        )

    parser.add_argument("execname", help="the name of the executable in cmssw project to run the jobs of")
    parser.add_argument("cfg",      help="the filename (relational path) of the cfg.py template for the jobs")
    parser.add_argument("dsets",    help="the filename (relational path) of the dsets json with dtag-dset targets for jobs")
    parser.add_argument("outdir",   help="the directory (relational path) for the jobs (FARM, cfg.py-s, input, output)")
    parser.add_argument("-d", "--dsets-dir",  help="use local directory for dsets files")
    parser.add_argument("--tausf",  help="turn on tau ID efficiency SF in cfg.py of jobs",
        action = "store_true")

    #if len(argv) != 5:
        #print("Usage:\njob_submit.py executable_filename cfg.py_template_filename dsets.json outdirname")
        #exit(1)

    # Process input

    #exec_name, cfg_templ_filename, dsets_json_filename, outdirname = argv[1:]

    args = parser.parse_args()

    #print(dir(args))
    #print(args.dsets_dir)
    #exit(1) # will crash

    exec_name = args.execname
    cfg_templ_filename = args.cfg
    dsets_json_filename = args.dsets
    outdirname = args.outdir
    outdirname = os.path.abspath(outdirname.strip()) + '/' # just in case

    # if retrieve dsets file info from local directory
    # or use DAS CLI client (failes from time to time)
    if args.dsets_dir:
        file_info_source = lambda dset: restore_dset_files_info(args.dsets_dir + '/', dset)
    else:
        file_info_source = get_dset

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

    #os.system('export X509_USER_PROXY=' + proxy_file)
    #my_env = os.environ.copy()
    #my_env["X509_USER_PROXY"] = proxy_file
    os.environ["X509_USER_PROXY"] = proxy_file

    print("Proxy is ready.")


    # finding local computing resourses
    hostname = commands.getstatusoutput("hostname -f")[1]
    print("hostname = %s" % hostname)
    hostname = '.'.join(hostname.split('.')[1:])
    print("hostname = %s" % hostname)
    # on lxplus hostnames are lxplusNNN.cern.ch
    #known_file_tiers = {"cern.ch": ("T2_CH_CERN", "root://eoscms//eos/cms/")}
    local_tier, local_file_server = known_file_tiers.get(hostname, ("CMS_DEFAULT_GLOBAL_XRD", "root://cms-xrd-global.cern.ch/"))

    n_files_per_job = 10

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
            #dset_files = get_dset_files(dset)
            #file_server = get_dset_site(dset)
            #dset_files, file_server = get_dset(dset)
            # file_info_source returns list of files of the dataset
            # and list of data tiers (file servers), which contain 100% files of this dataset
            dset_files, file_tiers = file_info_source(dset)
            if not dset_files:
                print("FAILED to fetch files on dset %s" % dset)
                print("dset files len  = %s" % len(dset_files))
                print("dset 100p tiers = %s" % file_tiers)
                continue

            # if there are no 100% full data tier, the default server is used
            if local_tier in file_tiers:
                file_server = local_file_server
            else:
                file_server = default_file_server
            # could do file_server = local_file_server if local_tier in file_tiers else default_file_server
            # but afraid for Python versions
            print("Chose file_server %s" % file_server)

            #dset_files = out_rows[3:]
            print("Found {} files. Splitting {} per job.".format(len(dset_files), n_files_per_job))


            #sites = sites_rows[3:]

            # the list of bsubs job submition commands
            dset_bsubs = []

            #make_jobs(dtag, dset_files, n_files_per_job)
            # make jobs per file chunks
            for i, job_chunk in enumerate(chunks(dset_files, n_files_per_job)):
                job_name = exec_name + '_' + dtag + '_' + str(i) + '_'

                # if the site of dataset is CERN -- can run directly on from eos
                #root://eoscms//eos/cms/
                # finding "physical" files via root://cms-xrd-global.cern.ch/
                input_files = ',\n'.join('"%s/%s"' % (file_server, s) for s in job_chunk)

                #if i < 5: print(input_files)

                job_cfg = cfg_templ.format(input = input_files, lumiMask = lumiMask, dtag = dtag, job_num = i, isMC = not isdata,
                                           outfile = outdirname + dtag + '_' + str(i) + '.root',
                                           outdir = outdirname + '/',
                                           project_dir = project_dir,
                                           withTauIDSFs = args.tausf)
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








