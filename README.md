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

The list of modules:

* `UserCode/llvv_fwk` from this repository: `https://github.com/cms2l2v/2l2v_fwk` -- it's used for POG definitions (Moriond17 etc), now lepton IDs and jet IDs and corrections are used in the ntuple output
* (The rest of the list of modules is comming.)

Overall log of actions:

	cd <CMSSW_X_X_X/src>
	cmsenv

	mkdir UserCode
	cd UserCode
	git clone git@github.com:xealits/ttbar-leptons-80X.git

	git clone https://github.com/cms2l2v/2l2v_fwk.git
	mv 2l2v_fwk llvv_fwk              # it was probably renamed at some moment

	# now remove some not used procedures in llvv_fwk, since they require more penedencies to compile
	# people use it to write analyses in the same repository with POG procedures
	# some of their analyses depend on things not needed here
	# at the moment (4 May 2017) I get error with
	# In file included from /.../CMSSW_8_0_25/src/UserCode/llvv_fwk/src/HiggsUtils.cc:1:0:
	# .../UserCode/llvv_fwk/interface/HiggsUtils.h:23:49: fatal error: ZZMatrixElement/MELA/interface/Mela.h: No such file or directory
	# the HiggsUtils are not used in ttbar-leptons-80X project (PatUtils, MacroUtils are used)
	# so just remove the .cc file
	mv llvv_fwk/src/HiggsUtils.cc llvv_fwk/src/HiggsUtils-cc-backup

	cd ttbar-leptons-80X
	scram b -j 9

	# back to src/ for more dependencies
	git cms-addpkg DataFormats/PatCandidates
	



CMSSW/test/arch to PATH for scripts
-----------------------------------

Also for convenience add your `<CMSSW>/test/<arch>` to `PATH`:
some small utilities are in `test/` directory and `scram` copies them to `test` of the CMSSW project,
but usual `cmsenv` doesn't add this directory to `PATH`.

Change the `cmssenv` from the usual `alias cmsenv='eval `scramv1 runtime -sh`'` to

	alias cmsenv='eval `scramv1 runtime -sh` && export PATH="$CMSSW_BASE/test/$SCRAM_ARCH:$PATH"'


