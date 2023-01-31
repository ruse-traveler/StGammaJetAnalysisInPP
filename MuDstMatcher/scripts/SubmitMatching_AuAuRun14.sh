#!/bin/bash
# 'SubmitMatching.sh'
# Derek Anderson
#
# Use this to create xml scripts and submit them via
# 'star-submit'.  This will automatically generate a
# list of files to analyze.  The variable 'Path' will
# is the directory where the desired files are.
#
# NOTE: requires the script 'GenerateXMLnew.sh' (to generate
#       the job description file)


# job parameters
#Path="/star/data105/embedding/AuAu_200_production_2014/Pythia6_pt9_11_100_20192901/P18ih.SL18h/2014"
Path="$1"  # DO NOT CHANGE
MuList="mudst.pt911dir$2.list"  # DO NOT REMOVE $2
GeList="geant.pt911dir$2.list"  # DO NOT REMOVE $2
Prefix="auau200r14pt911dir$2"  # DO NOT REMOVE $2
Label="matchWithMcInfo"
Suffix=".root"
Nevnt=1000
Nfile=1

# submission parameters
sim="false"
num="1"
cwd=$PWD
sub=$PWD"/submit"
ver="SL18h"  # DO NOT CHANGE
mac="bemcCalibMacro.C"  # DO NOT CHANGE
opd=$PWD"/output/"$Prefix
log=$PWD"/logs"
# do not modify
mu='\\\"$FILEBASENAME.MuDst.root\\\"'  # DO NOT CHANGE
ge='\\\"$filename\\\"'  # DO NOT CHANGE
lst=$GeList



# generate lists
printf "\n  Running submission script...\n"
declare -a Runs
declare -a Files

# modified for Run14 embedding sample
(( nRdir=0 ))
cd $Path
for dDir in {.*,*}; do
  if [ $dDir != "." ]; then
    if [ $dDir != ".." ]; then
      cd $dDir
      for rDir in {.*,*}; do
        if [ $rDir != "." ]; then
          if [ $rDir != ".." ]; then
            Runs[nRdir]=$(echo $dDir/$rDir)
            (( nRdir++ ))
          fi
        fi
      done  # end of run loop
      cd ..
    fi
  fi
done  # end of dir loop

# at most there will be 36 files in a run
(( nFiles=0 ))
for f in `seq 0 36`; do
  Files[nFiles]=$f;
  (( nFiles++ ))
done
printf "    Generated input lists...\n"



# loop over days
cd $cwd
touch $MuList
touch $GeList
(( nFiles=0 ))
# loop over runs
for run in ${Runs[@]}; do

  runDir=$Path"/"$run
  if [ ! -d $runDir ]; then
    continue
  fi

  # check if bad run
  if [[ $runDir == *"do_not_use"* ]]; then
    continue
  fi


  # loop over files
  # [01.30.2023] modified for Run14 emebdding sample 
  #for file in ${Files[@]}; do
  for file in $runDir/*.MuDst$Suffix; do

    #muFile=$Prefix"_"$run"_"$file".MuDst"$Suffix
    muFile=$file
    #muPath=$runDir"/"$muFile
    muPath=$file
    if [ ! -f $muPath ]; then
      continue
    fi

    #geFile=$Prefix"_"$run"_"$file".geant"$Suffix
    geFile=$file
    #gePath=$runDir"/"$geFile
    gePath=$file
    if [ ! -f $gePath ]; then
      continue
    fi

    # generate input list
    printf "$muPath" >> $MuList
    printf "$gePath" >> $GeList
    printf "\n" >> $MuList
    printf "\n" >> $GeList
    (( nFiles++ ))

  done  # end file loop
done  # end run loop
printf "    Generated file list '$GeList'...\n"
printf "      nFiles = $nFiles\n"


# generate output directory
./GenerateDir.sh $opd
printf "    Generated output directory '$opd'...\n"

# generate xml file
job=$Prefix"_"$Label
xml=$job".job.xml"
out='\\\"'$opd"/"'$FILEBASENAME.'$Label'.root\\\"'
arg=$(echo -e "\($Nevnt,$mu,$out,$Nfile,$ge\)")
./GenerateXML.sh $xml $sim $num $sub $ver $mac $arg $lst $log
printf "    Generated xml '$xml'...\n"

# generate log directory
if [ ! -d $log ]; then
  printf "  SubmitMatching.sh: Creating directory '$log'\n"
  mkdir $log
fi

# submit job
if [ ! -d $sub ]; then
  printf "  SubmitMatching.sh: Creating directory '$sub'\n"
  mkdir $sub
fi
cp $mac $sub
mv $xml $sub
mv $MuList $sub
mv $GeList $sub
cd $sub
star-submit $job".job.xml"
cd $cwd
printf "    Submited job '$job'...\n"


# delete input lists
unset Runs
unset Files

printf "  Finished submitting!\n\n"
