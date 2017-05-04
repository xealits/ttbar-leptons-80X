# very bad practice:
from job_handling import *
import yaml
import logging


# why __main__ if the whole file is script?
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = "Submit jobs to LSF cluster on lxplus.",
        epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
        )

    parser.add_argument("execname", help="the name of the executable in cmssw project to run the jobs of")
    parser.add_argument("defaults", help="the filename (relational path) of the YAML file with default configs")
    parser.add_argument("-c", "--configs", help="the filename (relational path) of the YAML file with additional configs")
    parser.add_argument("template",      help="the filename (relational path) of the cfg.py template for the jobs")
    parser.add_argument("-n", "--number-of-files", type=int, default=10, help="amount of files per job")
    parser.add_argument("dtags",    help="the filename (relational path) of the dtags yaml with dtag-dset targets for jobs (a dtag is a bunch of similar datasets, like _extN additions in MC or all data together)")
    parser.add_argument("-l", "--lumi-mask",    help="filename for the lumi-mask file to use for all jobs")
    parser.add_argument("outdir",   help="the directory (relational path) for the jobs (FARM, cfg.py-s, input, output)")
    parser.add_argument("--no-submit",  help="don't submit the generated jobs",
        action = "store_true")
    parser.add_argument("--test-configs",  help="don't generate and don't submit jobs, just check the config",
        action = "store_true")
    parser.add_argument("-d", "--dsets-dir",  help="use local directory for dsets files (small database of files in datasets and their current location)")
    parser.add_argument("--tausf",  help="turn on tau ID efficiency SF in cfg.py of jobs",
        action = "store_true")
    parser.add_argument("-s", "--file-server",  help='set file-server for the jobs, example: "root://eoscms//eos/cms/"')
    parser.add_argument("-m", "--more-configs", nargs='*', help='final configs of var=val form, overwriting the config files')

    #if len(argv) != 5:
        #print("Usage:\njob_submit.py executable_filename cfg.py_template_filename dsets.json outdirname")
        #exit(1)

    # Process input

    #exec_name, cfg_templ_filename, dsets_yaml_filename, outdirname = argv[1:]

    args = parser.parse_args()

    #print(dir(args))
    #print(args.dsets_dir)
    #exit(1) # will crash

    exec_name = args.execname
    dsets_yaml_filename = args.dtags
    outdirname = args.outdir
    outdirname = os.path.abspath(outdirname.strip()) + '/' # just in case

    # if retrieve dtags file info from local directory
    # or use DAS CLI client (failes from time to time)
    if args.dsets_dir:
        file_info_source = lambda dset: restore_dset_files_info(args.dsets_dir + '/', dset)
    else:
        file_info_source = get_dset

    with open(args.template) as t_f, open(dsets_yaml_filename) as d_f:
        #dsets = json.load(d_f)
        dtags = yaml.load(d_f)
        cfg_templ = t_f.read()

    with open(args.defaults) as cfg_defs:
        defaults  = yaml.load(cfg_defs)
    configs = defaults

    if args.configs:
        with open(args.configs) as cfg_defs:
            additional_configs  = yaml.load(cfg_defs)
        configs.update(additional_configs)

    if args.more_configs:
        # var=val list
        logging.info(args.more_configs)
        top_configs = dict([c.split('=') for c in args.more_configs])
        logging.info(top_configs)
        configs.update(top_configs)

    job_dir = outdirname + "/FARM/inputs/"
    job_outs = outdirname + "/FARM/outputs/"


    # Make job FARM directories

    logging.info("Making the task (job batch) area in {}".format(outdirname))

    try:
        os.makedirs(job_dir)
        os.makedirs(job_outs)
    except OSError, e:
        if e.errno == 17:
            logging.info("The directories exist, no need to make them.")
        else:
            raise

    project_dir = os.getcwd() + '/' # just in case
    logging.info("Project dir  is set to: %s" % project_dir)

    # could use list comprehension instead of this
    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        for i in xrange(0, len(l), n):
            yield l[i:i+n]

    # Create proxy for working with LSF jobs
    # Go through dtags.yaml
    # Create jobs.sh
    # Submit them

    if args.test_configs:
        logging.info("Config test: all configs are parsed, exiting.")
        exit(0)

    proxy_file = outdirname + "/FARM/inputs/x509_proxy"
    logging.info("Initializing proxy in " + proxy_file)
    #status, output = commands.getstatusoutput('voms-proxy-init --voms cms') # instead of voms-proxy-init --voms cms --noregen
    # os.system handles the hidden input well
    # subprocess.check_output(shell=True) and .Popen(shell=True)
    # had some issues
    status = os.system('voms-proxy-init --voms cms --out ' + proxy_file)

    if status !=0:
        logging.error("Proxy initialization failed: status = %s output =\n%s" % (str(status), output))
        exit(2)

    #os.system('export X509_USER_PROXY=' + proxy_file)
    #my_env = os.environ.copy()
    #my_env["X509_USER_PROXY"] = proxy_file
    os.environ["X509_USER_PROXY"] = proxy_file

    logging.info("Proxy is ready.")


    # finding local computing resourses
    hostname = commands.getstatusoutput("hostname -f")[1]
    logging.info("hostname = %s" % hostname)
    hostname = '.'.join(hostname.split('.')[1:])
    logging.info("hostname = %s" % hostname)
    # on lxplus hostnames are lxplusNNN.cern.ch
    #known_file_tiers = {"cern.ch": ("T2_CH_CERN", "root://eoscms//eos/cms/")}
    local_tier, local_file_server = known_file_tiers.get(hostname, ("CMS_DEFAULT_GLOBAL_XRD", "root://cms-xrd-global.cern.ch/"))

    n_files_per_job = args.number_of_files

    for dtag in dtags:
        dtag_name = dtag['dtag']
        logging.info("Submitting dtag " + dtag_name)

        isdata = dtag.get('isMC', False)
        if args.lumi_mask:
            lumiMask = args.lumi_mask
        else:
            lumiMask = dtag.get('lumiMask', '')

        dtag_files = [] # the physicl paths for the files of the dtag jobs
        for dset in dtag['dsets']:
            logging.info(dset)

            # get files of the dset
            # and list of data tiers (file servers), which contain 100% files of this dataset
            # combine them into list of "physical paths" and store all lists together for the whole dtag
            dset_files, file_tiers = file_info_source(dset)
            if not dset_files:
                logging.error("FAILED to fetch files on dset %s\n dset files len  = %s\n dset 100p tiers = %s" % (dset, len(dset_files), file_tiers))
                continue

            if args.file_server:
                file_server = args.file_server
            # if there are no 100% full data tier, the default server is used
            elif local_tier in file_tiers:
                file_server = local_file_server
            else:
                file_server = default_file_server
            # could do file_server = local_file_server if local_tier in file_tiers else default_file_server
            # but afraid for Python versions
            logging.info("Chose file_server %s" % file_server)

            # now we have the list of filenames
            # and the best guess for the server where to find them (either EOS of cern of cmsglobal)
            # the input file URL is "server/file"
            for f in dset_files:
                dtag_files.append("%s/%s" % (file_server, f))

        logging.info("Found {} files. Splitting {} per job.".format(len(dtag_files), n_files_per_job))

        # the list of bsubs job submition commands
        dtag_bsubs = []

        #make_jobs(dtag, dset_files, n_files_per_job)
        # make jobs per file chunks
        for i, job_chunk in enumerate(chunks(dtag_files, n_files_per_job)):
            job_name = exec_name + '_' + dtag_name + '_' + str(i) + '_'

            # if the site of dataset is CERN -- can run directly on from eos
            #root://eoscms//eos/cms/
            # finding "physical" files via root://cms-xrd-global.cern.ch/
            input_files = ',\n'.join(f for f in job_chunk)

            #if i < 5: print(input_files)

            # creating the _cfg.py file (which is just a container of parameters for the job)
            # first configs are filled with defaults and for each jobs
            # they are filled with parameters passed on command line
            configs.update({'input': input_files, 'lumiMask': lumiMask, 'dtag': dtag_name, 'job_num': i, 'isMC': not isdata,
                            'outfile'      : outdirname + '/' + dtag_name + '_' + str(i) + '.root',
                            'outdir'       : outdirname + '/',
                            'project_dir'  : project_dir,
                            'withTauIDSFs' : args.tausf})
            job_cfg = cfg_templ.format(**configs)

            job_cfg_filename = outdirname + '/' + dtag_name + '_' + str(i) + '_cfg.py'
            with open(job_cfg_filename, 'w') as cfg_file:
                cfg_file.write(job_cfg)

            # creating the job.sh file with is then submitted to the batch system
            job_sh = job_template.format(project_dir = project_dir, exec_name = exec_name, job_cfg = job_cfg_filename)
            job_sh_filename = outdirname + '/FARM/inputs/' + dtag_name + '_' + str(i) + '.sh'
            with open(job_sh_filename, 'w') as sh_file:
                sh_file.write(job_sh)
                os.system("chmod 777 " + sh_file.name)

            # record the submition command for this job (all these commands are written to a bsub.sh file, which is sourced (executed) to resubmit the jobs)
            job_stdout = outdirname + '/FARM/outputs/' + dtag_name + '_' + str(i) + '.cout'
            dtag_bsubs.append(job_command_template.format(job_name = job_name, job_stdout = job_stdout, jobsh = job_sh_filename))

        logging.info("Created {} jobs. Ready to submit.".format(len(dtag_bsubs)))

        # write the submition commands to the shell script file
        with open(outdirname + '/FARM/inputs/' + dtag_name + '_' + 'bsub.sh', 'w') as bsub_file:
            bsub_file.write('\n'.join(dtag_bsubs))
            logging.info("Wrote the bsub commands to {}".format(bsub_file.name))

        logging.info("Ready to submit.")
        #print("No submit in the test run.")
        #print("printing the bsubs instead")

        # if set to submit the jobs right away then do so
        if args.no_submit:
            logging.info("NO SUBMIT, as requested on command line.")
        else:
            for bsub in dtag_bsubs:
                #print(bsub)
                logging.info(commands.getstatusoutput(bsub))

