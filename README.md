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
* (it seems the rest of dependencies are CMSSW built-ins...)

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

	# and done! (it even compiled on CMSSW_8_0_26, which has "latest greatest MET filters" and I need to switch to them..)
	# I'm testing if anything else is really needed
	# like this stuff:
	# back to src/ for more dependencies
	#git cms-addpkg DataFormats/PatCandidates
	

alternatively you can run ./todo.pl which will install all necessary codes


CMSSW/test/arch to PATH for scripts
-----------------------------------

Also for convenience add your `<CMSSW>/test/<arch>` to `PATH`:
some small utilities are in `test/` directory and `scram` copies them to `test` of the CMSSW project,
but usual `cmsenv` doesn't add this directory to `PATH`.

Change the `cmssenv` from the usual `alias cmsenv='eval `scramv1 runtime -sh`'` to

	alias cmsenv='eval `scramv1 runtime -sh` && export PATH="$CMSSW_BASE/test/$SCRAM_ARCH:$PATH"'






Submition of jobs, working with data
=======================================

Datasets (called "dset" throughout the jobing code) are treated in bunches called "dtag".
A "dtag" contains masically the same datasets
-- mainly "-ext" Monte-Carlo datasets, but also one might join several Data runs in 1 dtag.
Now 1-to-1 pairs of dtag-dset are used, but the -ext MC is to come.

Jobs are named "dtag_<number of job>" and produce corresponding "<job>.<subchannel>.root", "<job>.json" (luminosities), "<job>.job_done" (for control of resubmit etc).
    
Thus datasets in a dtag should have the same cross-sections for correct further processing.

There is a small local "database" made with `store_dsets_info.py` which runs `das_client` and stores locally files of given datasets and their current location.
So that one doesn't need to pull this list of files on every submition of the jobs.
The "location" is the tier on which 100% of dataset files (not blocks) are present.
Currently only CERN's EOS (`T2_CERN_CH`) is distinguished by job submition utility (`job_submit.py`).
But it is easily extendable.





# NTuple production

The NTuple interface is in `interface/ntupleOutput.h`.
It is defined with list of macro, which are unpacked by preprocessor in some local namespace of your program
and attached to a TTree you've provided.
The comments in the `ntupleOutput.h` describe it better.
Roughly, you use it like this in your `main`:

	int main()
	{
	...
	TTree NT_output_ttree("reduced_ttree", "TTree with reduced event data");
	// to create branches in the new TTree:
	#define NTUPLE_INTERFACE_CREATE
	// for attaching to existing branches define NTUPLE_INTERFACE_OPEN
	// unpack the interface
	#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput.h"

	// the code of the interface is unpacked
	// now there are all the NT_branchName variable in the name space
	// and their pointers are passed to branches

	for (in the loop of events)
		{
		RESET_NTUPLE_PARAMETERS // macro defaults all parameters
		NT_tauID = event[i].taus[0].pdgId();
		NT_output_ttree.Fille();
		}
	}

To add new objects to the interface use this macro:

	#define OBJECT_in_NTuple(NTuple, Class, Name)   Class   NT_##Name; NTuple.Branch(#Name, #Class, &NT_##Name)

as in following:

	OBJECT_in_NTuple(OUTNTUPLE, Pat::Tau, tau1);

-- but ROOT needs to "know" the class to be able to add as Branch.
I have not tested yet which classes pass by default.
They say to add a class to root you run this:

	gROOT->ProcessLine("#include <vector>")

-- and now you can save vectors of floats/ints etc.

Also, if you add a new object, consider adding it to `RESET_NTUPLE_PARAMETERS`.
You need to add a command which resets/clears the object for each iteration over events:

	NT_myBranchNameA.clear(); \
	NT_myBranchNameB.reset();

-- unless you don't need to clear the object on for each event.
(Uncleared value might get stored into the TTree.)



MISC
====

Just to keep it around:

    http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_8_0_25/doc/html/d1/de9/classpat_1_1Tau.html

