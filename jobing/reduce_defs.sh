
# No copy, access the scripts via path in env var, obtained by readlink -e

SCRIPTS_SUM_COUNTERS=`readlink -e job-reduce/processing-counters_merging-sets.R`
SCRIPTS_SUM_DISTR=`readlink -e job-reduce/processing-distributions_merging-sets.R`
SCRIPTS_STITCH_WEIGHTFLOW=`readlink -e job-reduce/processing_stitch_weightflow.R`

# No copy, access output names via env var:
# store full path to the job-reduce/ dir and access files in there

JOBREDUCE_DIR=`readlink -e job-reduce/`



# +++++++++ merge procedures

function Merge_distr {
   # Description: 
   # 2 input parameters
   #   $1 distr name
   #   $2 distr header name
   distrs_headername=$2

   { # adding the red collor to error stream
   echo Merge_distr $1 $2
   echo "type,dtag,job_num," `grep -m 1 --no-filename "^header,$distrs_headername," ../*.csv | head -n 1 | sed "s/^header,$distrs_headername,//"` > $1
   grep --no-filename "^$1:content," ../*csv | sed "s/^$1:content,//" >> $1 
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done) && echo "# Merged distr >"
}  


function Merge_counter {
   # Description: 
   # 1 input parameter
   #   $1 counter name

   { # adding the red collor to error stream
   echo Merge_counter $1
   echo "type,dtag,job_num,$1" > $1
   grep --no-filename "^$1," ../*csv | sed "s/^$1,//" >> $1 
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done) && echo "# Merged counter >"
}





# +++++++++ sum procedures
function Sum_distr {
   # Description:
   echo Sum_distr $1
   input=../"$1"
   output="$1"

   { # adding the red collor to error stream
   # ./processing-distributions_merging-sets.R $input $output
   $SCRIPTS_SUM_DISTR $input $output
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done) && echo "# Summed distr >"
}


function Sum_counter {
   # Description:
   echo Sum_counter $1
   input=../"$1"
   output="$1"

   { # adding the red collor to error stream
   #./processing-counters_merging-sets.R $input $output
   $SCRIPTS_SUM_COUNTERS $input $output
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done) && echo "# Summed counter >"
}



# ============ both procedures

function Reduce_distr {
   # Description: 
   # 3 input parameters
   #   $1 counter name
   #   $2 distr header name
   #   $3 directory with jobs outputs

   { # adding the red collor to error stream
   echo in $3 reducing $1
   cd $3/merged-sets
   Merge_distr $1 $2
   cd jobsums
   Sum_distr $1
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done)
}



function Reduce_counter {
   # Description: 
   # 2 input parameters
   #   $1 counter name
   #   $2 directory with jobs outputs

   { # adding the red collor to error stream
   echo in $2 reducing $1
   cd $2/merged-sets
   Merge_counter $1
   cd jobsums
   Sum_counter $1
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done)
}


