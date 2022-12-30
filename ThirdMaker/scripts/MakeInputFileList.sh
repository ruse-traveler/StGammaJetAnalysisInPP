#!/bin/bash
# 'MakeInputFileList.sh'
#
# Need a list of files? Use this.

list="files.all2.d20m3y2017.list"

get_file_list.pl -keys node,path,filename -cond filename~st_gamma,trgsetupname=production2009_200Gev_Hi,production=P11id,filetype=daq_reco_MuDst,filename~st_gamma,storage!=HPSS -limit 0 > $list

