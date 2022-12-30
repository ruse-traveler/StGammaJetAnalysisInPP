#!/bin/bash
# 'MakeTreeList.sh'
#
# Use this to make lists of 'StThirdMaker' and
# 'StThirdJetMaker' output.

tuIn="/star/data01/pwg/dmawxc/JetReco_pp/FullJetTreeTest/tuple.4*"
trIn="/star/data01/pwg/dmawxc/JetReco_pp/FullJetTreeTest/NewBuild/tree2.25.*"
tuOut="files.tuple.d3m3y2017.list"
trOut="files.tree2.25.d3m3y2017.list"

ls $tuIn > $tuOut
ls $trIn > $trOut
