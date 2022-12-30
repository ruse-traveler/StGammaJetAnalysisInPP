#! /bin/csh
# 'TriggerResolutionStudy_RunStarsim.csh'
# Derek Anderson
# 08.30.2020
#
# Runs Particle Gun framework:
#   1st) simulates particles via starsim
#        and produces .fzd file
#   2nd) generates mudst using .fzd as
#        input
# Adapted from code from Robert Licenik
# and Maria Zurek.


# set events
if ( $#argv == 0 ) then
  set nEvents=10
else
  set nEvents=$argv[1]
endif

# announce start
echo "====== Begin Starsim ====="
echo "  Setting nEvents: $nEvents"


# set random seed
set seed=`date +%N`
echo "  Setting seed to: $seed"

# set output name
set outputFile="pionGridTest_onlyOneEtaSite_pin8pic2"
echo "  Setting output name: $outputFile"

# remove old files
rm -f kinematics_$outputFile.starsim.fzd


# set star library for starsim
set starSimLib=SL20a
echo "  Settting star library: $starSimLib"
starver $starSimLib

# run starsim
echo "  Running starsim for $nEvents events"
root4star -l -b -q 'TriggerResolutionStudy_RunStarsim.C('$nEvents', '$seed', "'$outputFile'")'


# set star library for bfc
set starBfcLib=SL11d_embed
echo "  Setting star library: $starBfcLib"
starver $starBfcLib

# run bfc
echo "  Running bfc reconstruction for $nEvents events"
root4star -l -b -q 'TriggerResolutionStudy_RunBfc.C('$nEvents', "kinematics_'$outputFile'.starsim.fzd")'

# announce end
echo "====== End  Starsim ======"

# End -------------------------------------------------------------------------
