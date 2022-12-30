#!/bin/bash
# 'resubmit.sh'
#
# Use this to re-submit jobs.

xml=87311E8F310B172664D4DAA9A85E46BC.session.xml
list1=348-351
list2=848-866
list3=3336,3337

star-submit -r $list1 $xml
star-submit -r $list2 $xml
star-submit -r $list3 $xml
