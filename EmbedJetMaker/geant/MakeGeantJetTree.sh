#!/bin/bash
# 'MakeGeantJetTree.sh'
# Derek Anderson
# 10.03.2017
#
# Use this to run 'MakeGeantJetTree.C'
# in batch mode.

root -b -q MakeGeantJetTree.C\("true"\)
