#! /usr/bin/perl
use Cwd;
use POSIX;
use POSIX qw(strftime);

#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];
$UserID= POSIX::cuserid();
$UserIDCern=$UserID;
$UserDir="";
$CMSSWRel = "8_0_25";
$ARCH='slc6_amd64_gcc530';
$GitSource="xealits"


if($ARGV[0] eq "--help" || $ARGV[0] eq ""){
    printf("\nThe installation follows the recomendation given in the description");
    printf("\nof  for the latest CMSSW release: https://github.com/LLRCMS/LLRHiggsTauTau");
    printf("\n\nThis code requires one input option. The syntax is: ./todo.pl [OPTION]");
    printf("\nPlease choose from the following options:\n");
    printf("\n./todo.pl --help                                   Prints this message\n");
    printf("\n./todo.pl --TopNtuple <TopNtupleDir>               Setups up TopNtuple and gives instructions for submission");
    printf("\n                                                   <TauNtupleDir> location of CMSSW for TauNtuple.");
    printf("\n                                                      Optional Commmads: ");
    printf("\n                                                    --CMSSWRel <CMSSW release #> The CMSSW release you want to use Default: $CMSSWRel");
    printf("\n                                                    --GIT git user to take  ttbar-leptons-80X repo (if forked). Default: xealits ( master  branch) \n\n\n");

    exit(0);  
}
 my $dir = getcwd;
$time= strftime("%h_%d_%Y",localtime);

for($l=1;$l<$numArgs; $l++){
    if($ARGV[$l] eq "--GIT"){
	$GitSource=$ARGV[$l+1];
    }
} 
if( $ARGV[0] eq "--TopNtuple"){
    $currentdir=getcwd;
    if($ARGV[1] ne ""){
	$basedir=$ARGV[1];
    }
    else{
	printf("\nWorkingDir for CMSSW is required. Please follow the syntax:./todo.pl --TopNtuple <TopNtupleDir> ");
	printf("\nFor more details use: ./todo --help\n"); 
	exit(0);
    }
    printf("\nWorkingDir for CMSSW: $basedir");
    printf("\nCurrentDir is: $currentdir \n");
    system(sprintf("rm Install_TopNtuple_$time"));
    system(sprintf("echo \"cernlib-use root\" >> Install_TopNtuple_$time"));
    system(sprintf("echo \"export SCRAM_ARCH=\\\"$ARCH\\\"\" >> Install_TopNtuple_$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Install_TopNtuple_$time"));
     system(sprintf("echo \"mkdir $basedir\" >>  Install_TopNtuple_$time")); 
    system(sprintf("echo \"cd $basedir\" >>  Install_TopNtuple_$time")); 
    system(sprintf("echo \"cmsrel CMSSW_$CMSSWRel\" >>  Install_TopNtuple_$time")); 
    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \"cmsenv;\" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \"git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit;    \" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \"cd HiggsAnalysis/CombinedLimit ; git checkout 81x-root606-integration;  \" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \" source env_standalone.sh ;   \" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \"  make -j 8; make;  \" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \"  cd ../../;  \" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \" mkdir UserCode; cd UserCode;\" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \" git clone https://github.com/$GitSource/ttbar-leptons-80X.git;\" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \" git clone https://github.com/cms2l2v/2l2v_fwk.git;\" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \"mv 2l2v_fwk llvv_fwk    ;\" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \" mv llvv_fwk/src/HiggsUtils.cc llvv_fwk/src/HiggsUtils-cc-backup  ;\" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \" cd ttbar-leptons-80X  ;\" >> Install_TopNtuple_$time")); 
    system(sprintf("echo \"  scram b -j 9 ;\" >> Install_TopNtuple_$time")); 



    # print Instructions
    printf("\n\nInstructions:");
    printf("\n source   Install_TopNtuple_$time  \n\n");
 

}
 
