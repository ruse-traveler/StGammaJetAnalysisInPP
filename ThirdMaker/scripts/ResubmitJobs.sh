#!/bin/bash
# 'ResubmitJobs.sh'
#
# Use this to resubmt some jobs...

jobz="all"
sess="C27C92D9FEE17EE88FCFF22D7ECE052B.session.xml"

star-submit -kr $jobz $sess
