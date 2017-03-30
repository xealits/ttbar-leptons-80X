TTbar leptons in CMSSW80X
=========================

The analysis/event selection of ttbar-to-leptons decay in CMSSW80X.





Installation
============

Just check out the repository in your CMSSW project:

	cd <CMSSW_X_X_X/src>
	mkdir UserCode
	cd UserCode
	git clone git@github.com:xealits/ttbar-leptons-80X.git

And add all used CMSSW modules in `<CMSSW_X_X_X/src>`.
(The list of modules is comming.)

Also for convenience add your `<CMSSW>/test/<arch>` to `PATH`:
some small utilities are in `test/` directory and `scram` copies them to `test` of the CMSSW project,
but usual `cmsenv` doesn't add this directory to `PATH`.

Change the `cmssenv` from the usual `alias cmsenv='eval `scramv1 runtime -sh`'` to

	alias cmsenv='eval `scramv1 runtime -sh` && export PATH="$CMSSW_BASE/test/$SCRAM_ARCH:$PATH"'


