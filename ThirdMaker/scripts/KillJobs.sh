#!/bin/bash
# 'KillJobs.sh'
#
# Use this to kill jobs.

xml=C27C92D9FEE17EE88FCFF22D7ECE052B.session.xml
list1=all
#list2=848-866
#list3=3336,3337

star-submit -k $list1 $xml
#star-submit -k $list2 $xml
#star-submit -k $list3 $xml
