# 'PushToHpss.sh'
#
# Use this to copy data to HPSS.

date="10Oct2016"

echo "Copying to HPSS! Today is $date"

htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia01.jets."$date".tar ./Pythia01* >& Pythia01.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia02.jets."$date".tar ./Pythia02* >& Pythia02.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia05.jets."$date".tar ./Pythia05* >& Pythia05.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia07.jets."$date".tar ./Pythia07* >& Pythia07.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia08.jets."$date".tar ./Pythia08* >& Pythia08.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia09.jets."$date".tar ./Pythia09* >& Pythia09.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia10.jets."$date".tar ./Pythia10* >& Pythia10.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia11.jets."$date".tar ./Pythia11* >& Pythia11.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia15.jets."$date".tar ./Pythia15* >& Pythia15.hpss.log & 
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia18.jets."$date".tar ./Pythia18* >& Pythia18.hpss.log & 
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia19.jets."$date".tar ./Pythia19* >& Pythia19.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia20.jets."$date".tar ./Pythia20* >& Pythia20.hpss.log & 
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia21.jets."$date".tar ./Pythia21* >& Pythia21.hpss.log & 
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/Pythia22.jets."$date".tar ./Pythia22* >& Pythia22.hpss.log & 
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetReconstruction/tests.jets."$date".tar tests/* >& tests.hpss.log & 
