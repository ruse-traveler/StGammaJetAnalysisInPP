#!/bin/bash
# 'TriggerResolutionStudy_Submit.sh'
# Derek Anderson
# 09.04.2020
#
# Submits jobs to run Particle Gun
# framework. Adapted from code from
# Robert Licenik and Maria Zurek.
#
# Note: submitdir=`/bin/pwd -P`

# set parameters
starversion=SL11d_embed
submitdir=$PWD
echo $submitdir
prodId=`date +%F_%H-%M`
submitoutdir=submitout
simoutdir=simdata
recooutdir=recodata

# no. of events and random seed
nevents=100
random=$RANDOM

# additional parameters (uncomment as necessary)
#simulator=("pythia8" "pythia8tune" "pythia6" "pythia6perugia0")
#simulator=("pythia8" "pythia8tune")
#simulator=("pythia6" "pythia6perugia0")
#evetype=("minbias" "promptjpsi")
#evetype=("minbias" "sngd0") #data2, old kine cuts for MTD and 5 D0 per event
#evetype=("sngd0") #data3
#evetype=("minbias" "sngd0") #data4,5


# create output directories
mkdir -p ./SubmitInfo/

  echo "output directories"
  echo ./$submitoutdir/${prodId}
  echo ./$simoutdir/${prodId}
  echo ./$recooutdir/${prodId}

  mkdir -p ./$submitoutdir/${prodId}
  mkdir -p ./$simoutdir/${prodId}
  mkdir -p ./$recooutdir/${prodId}


# submit jobs
star-submit-template -template TriggerResolutionStudy_Submit.xml -entities submitdir=$submitdir,submitoutdir=$submitoutdir,simoutdir=$simoutdir,recooutdir=$recooutdir,nevents=$nevents,random=$random,prodId=$prodId,starversion=$starversion

# End -------------------------------------------------------------------------
