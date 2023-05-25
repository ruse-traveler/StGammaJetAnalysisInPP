# 'ExecuteAnalysis.sh'
# Derek Anderson
#
# This is for performing the gamma+jet analysis in batch mode

echo "  Executing analysis..."
nohup ./GammaJetExecute.exe GammaJetExecute.cmnd >& GammaJetExecute.out &
echo "  Analysis executed!"
