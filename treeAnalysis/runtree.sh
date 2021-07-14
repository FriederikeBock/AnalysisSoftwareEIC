input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB")'
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_pTHard5")'
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB_pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB_pTHard5")'


input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-TREXTOUT_epMB.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-TREXTOUT_epMB")'
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-TREXTOUT_pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-TREXTOUT_pTHard5")'
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-TREXTOUT_epMB_pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-TREXTOUT_epMB_pTHard5")'


if [ $1 = "DRpion" ]; then
    input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/tree_CALOSTANDALONE_SimplePhoton.root
    root -x -l -b -q 'treeProcessing.C("'$input'","CALOSTANDALONE_SimplePhoton",true,true,true,true,true,true,true,true,true,true,true,true,false, true, true, -1, 0, 1, false)'
fi

if [ $1 = "calogamma" ]; then
    input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/tree_CALOSTANDALONE_SimplePhoton.root
    root -x -l -b -q 'treeProcessing.C("'$input'","CALOSTANDALONE_SimplePhoton",true,true,true,true,true,true,true,true,true,true,true,true,false, true, true, -1, 0, 1, false)'
elif [ $1 = "calopion" ]; then
    input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/tree_CALOSTANDALONE_SimplePion.root
    root -x -l -b -q 'treeProcessing.C("'$input'","CALOSTANDALONE_SimplePion",true,true,true,true,true,true,true,true,true,true,true,true,false, true, true, -1, 0, 1, false)'
elif [ $1 = "fhcalpion" ]; then
    input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/tree_FHCALSTANDALONE_SimplePion.root
    root -x -l -b -q 'treeProcessing.C("'$input'","FHCALSTANDALONE_SimplePion",true,true,true,true,true,true,true,true,true,true,true,true,false, true, true, -1, 0, 1, false)'
fi

# MODULAR
if [ $1 = "1" ]; then
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB_pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","MODULAR_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true)'
elif [ $1 = "2" ]; then
input=/media/nschmidt/external3/EIC_outputs2/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250MBandpTh.root
root -x -l -b -q 'treeProcessing.C("'$input'","MODULAR_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true)'
elif [ $1 = "3" ]; then
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-TREXTOUT_epMB_pTHard5.root
root -x -l -b -q 'treeProcessing.C("'$input'","MODULAR_ALLSILICON-TREXTOUT_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true)'
elif [ $1 = "4" ]; then
input=/media/nschmidt/external3/EIC_outputs2/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT-HC2xEC2x_e10p250pTHard5.root
root -x -l -b -q 'treeProcessing.C("'$input'","MODULAR_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT-HC2xEC2x_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true, true, true, -1, 0, 2)'

# BEAST MODULAR
elif [ $1 = "5" ]; then
input=/media/nschmidt/external3/EIC_outputs2/BeastModular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_e10p250pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","BEAST_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true)'
elif [ $1 = "6" ]; then
input=/media/nschmidt/external3/EIC_outputs2/BeastModular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250MBandpTh.root
root -x -l -b -q 'treeProcessing.C("'$input'","BEAST_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250MBandpTh",true,true,true,true,true,true,true,true,true,true,true,true,true)'
elif [ $1 = "7" ]; then
input=/media/nschmidt/external3/EIC_outputs2/BeastModular/output_ALLSILICON-TREXTOUT_e10p250pTHard5.root
root -x -l -b -q 'treeProcessing.C("'$input'","BEAST_ALLSILICON-TREXTOUT_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true)'
fi

if [ $1 = "newinput_inlay" ]; then
    input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-ETTL-CTTL-INNERTRACKING-DRCALO-FwdSquare-FTTLDRC-DRTungsten-GEOMETRYTREE-TRACKEVALHITS_MBMC/G4EICDetector_eventtree.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'geometry.root'","NEWINPUT_GEOMETRY",true,true,false,  false,true,true,-1,0,false,0)'
fi

if [ $1 = "newinput_pythia" ]; then
    #input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-ETTL-CTTL-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-ASYM-FTTLS3LC-GEOMETRYTREE-TEST/G4EICDetector_eventtree.root
    #geometry=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-ETTL-CTTL-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-ASYM-FTTLS3LC-GEOMETRYTREE-TEST/geometry.root
    input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-ASYM-GEOMETRYTREE-PYTHIA-1x/G4EICDetector_eventtree.root
    geometry=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-ASYM-GEOMETRYTREE-PYTHIA-1x/geometry.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'$geometry'","NEWINPUT_TMSTUD_PYTHIA",true,false  ,true,true,-1,0,false,0)'
fi


if [ $1 = "newdet_pythia" ]; then
    input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-LFHCAL-GEOMETRYTREE-EEMCH-BECAL-PYTHIA-1/G4EICDetector_eventtree.root
    geometry=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-LFHCAL-GEOMETRYTREE-EEMCH-BECAL-PYTHIA-1/geometry.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'$geometry'","NEWDETS_PYTHIA",true,false  ,true,true,-1,0,false,0)'
fi


if [ $1 = "newdet_pythia_fix" ]; then
    input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-LFHCAL-GEOMETRYTREE-EEMCH-BECAL-PYTHIA-EVALEIC-lll-tst2/G4EICDetector_eventtree.root
    geometry=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-LFHCAL-GEOMETRYTREE-EEMCH-BECAL-PYTHIA-EVALEIC-lll-tst2/geometry.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'$geometry'","PYTHIA_EVALPROJFIX",true,false  ,true,true,-1,0,false,0)'
fi


if [ $1 = "newinput_singlep" ]; then
    input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-ASYM-GEOMETRYTREE-SINGLEPIELEC/G4EICDetector_eventtree.root
    geometry=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-ASYM-GEOMETRYTREE-SINGLEPIELEC/geometry.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'$geometry'","NEWINPUT_TMSTUD_SINGLEP",true,false  ,true,true,-1,0,false,0)'
fi

