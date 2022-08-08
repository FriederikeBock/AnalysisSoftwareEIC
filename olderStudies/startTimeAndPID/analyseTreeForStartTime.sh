#!/bin/bash
date

detSetup="ALLSILICON-FTTLSE2LC-ETTL-CTTLSE1-ACLGAD-TREXTOUT_epMB"
root -l -b -q 'analyseTreeForStartTime.C+("'$detSetup'", "", "'dataList/$detSetup.list'", 0, 1, 1)'

#detSetup="ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB"
#root -l -b -q 'analyseTreeForStartTime.C+("'$detSetup'", "", "'dataList/$detSetup.list'", 0, 1, 1)'

#detSetup="ALLSILICON-FTTLSE2LC-ETTL-CTTLSE1-ACLGAD-TREXTOUT_epMB"
#root -l -b -q 'analyseTreeForStartTime.C+("'$detSetup'", "", "'dataList/$detSetup.list'", 0, 1, 0)'
