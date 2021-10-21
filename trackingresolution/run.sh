#root -x -l -q -b trackingreso.C
#root resolutionJETS.C
if [ $1 == "Nico" ]; then
basepath=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing_centralProduction_raymond
trkrid=1

#root -x -l -q -b 'trackingeffi.C("'$basepath'/production-singleElectron-p-0-to-20/output_TRKEFF.root","pdf","elec-centprod",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/production-singlePion-p-0-to-20/output_TRKEFF.root","pdf","pion-centprod",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption5/output_TRKEFF.root","pdf","elec-geo5",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption6/output_TRKEFF.root","pdf","elec-geo6",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption5/output_TRKEFF.root","pdf","pion-geo5",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption6/output_TRKEFF.root","pdf","pion-geo6",'$trkrid')'

root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-pythia8-10x100-q2-1-to-100-geoOption5/output_TRKEFF.root","pdf","pythia-centprod-geo5-q2-1-100",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-pythia8-10x100-q2-1-to-100-geoOption6/output_TRKEFF.root","pdf","pythia-centprod-geo6-q2-1-100",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-pythia8-10x100-q2-100-geoOption5/output_TRKEFF.root","pdf","pythia-centprod-geo5-q2-100",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-pythia8-10x100-q2-100-geoOption6/output_TRKEFF.root","pdf","pythia-centprod-geo6-q2-100",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption5/output_TRKEFF.root","pdf","elec-geo5",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption6/output_TRKEFF.root","pdf","elec-geo6",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption5/output_TRKEFF.root","pdf","pion-geo5",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption6/output_TRKEFF.root","pdf","pion-geo6",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/production-singleElectron-p-0-to-20/output_TRKEFF.root","pdf","elec-centprod",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/production-singlePion-p-0-to-20/output_TRKEFF.root","pdf","pion-centprod",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/production-pythia8-10x100-q2-1-to-100/output_TRKEFF.root","pdf","pythia-centprod-q2-1-100",'$trkrid')'
root -x -l -q -b 'trackingeffi.C("'$basepath'/production-pythia8-10x100-q2-100/output_TRKEFF.root","pdf","pythia-centprod-q2-100",'$trkrid')'

# trkrid=2
# trksource=0
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-pythia8-10x100-q2-1-to-100-geoOption5/output_TRKRS.root","pdf","pythia-centprod-geo5-q2-1-100", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-pythia8-10x100-q2-1-to-100-geoOption6/output_TRKRS.root","pdf","pythia-centprod-geo6-q2-1-100", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-pythia8-10x100-q2-100-geoOption5/output_TRKRS.root","pdf","pythia-centprod-geo5-q2-100", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-pythia8-10x100-q2-100-geoOption6/output_TRKRS.root","pdf","pythia-centprod-geo6-q2-100", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption5/output_TRKRS.root","pdf","elec-geo5", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption6/output_TRKRS.root","pdf","elec-geo6", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption5/output_TRKRS.root","pdf","pion-geo5", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption6/output_TRKRS.root","pdf","pion-geo6", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/production-singleElectron-p-0-to-20/output_TRKRS.root","pdf","elec-centprod", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/production-singlePion-p-0-to-20/output_TRKRS.root","pdf","pion-centprod", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/production-pythia8-10x100-q2-1-to-100/output_TRKRS.root","pdf","pythia-centprod-q2-1-100", true,'$trkrid','$trksource')'
# root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/production-pythia8-10x100-q2-100/output_TRKRS.root","pdf","pythia-centprod-q2-100", true,'$trkrid','$trksource')'

fi 


#make effi with other trackset


if [ $1 == "Fredi" ]; then
  if [ $2 == "tree" ]; then
    cd ../treeAnalysis/
    geometryfile=/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root
#     root -b -x -l -q 'treeProcessing.C("LYSO-1TTL-muRwellRes-pythia8Q2.txt","'$geometryfile'","LYSO-1TTL-muRwellRes-pythia8Q2-TRonly",false,false,true,false,-1,0,false,0,"anti-kt",0.5,30,false)'
#     root -b -x -l -q 'treeProcessing.C("LYSO-1TTL-muRwellRes-pythia8.txt","'$geometryfile'","LYSO-1TTL-muRwellRes-pythia8-TRonly",false,false,true,false,-1,0,false,0,"anti-kt",0.5,30,false)'
#     root -b -x -l -q 'treeProcessing.C("LYSO-1TTL-nomRes-pythia8Q2.txt","'$geometryfile'","LYSO-1TTL-nomRes-pythia8Q2-TRonly",false,false,true,false,-1,0,false,0,"anti-kt",0.5,30,false)'
#     root -b -x -l -q 'treeProcessing.C("LYSO-1TTL-nomRes-pythia8.txt","'$geometryfile'","LYSO-1TTL-nomRes-pythia8-TRonly",false,false,true,false,-1,0,false,0,"anti-kt",0.5,30,false)'
    root -b -x -l -q 'treeProcessing.C("TTL-def_HFJetsQ100.txt","'$geometryfile'","2TTL-nomRes-pythia8Q100-TRonly",-1,false,false,false,false,0,0,"anti-kt",0.5,30,false)'
    root -b -x -l -q 'treeProcessing.C("TTL-def_HFJetsQ1-100.txt","'$geometryfile'","2TTL-nomRes-pythia8Q1-100-TRonly",-1,false,false,false,false,0,0,"anti-kt",0.5,30,false)'
    
    cd ../trackingresolution/
  fi
  if [ $2 == "extract" ]; then
    basepath=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/LYSO-1TTL-muRwellRes-pythia8Q2-TRonly/output_TRKRS.root","pdf","1TTLmuRWellResMinQ2_10",false,0,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/LYSO-1TTL-muRwellRes-pythia8Q2-TRonly/output_TRKRS.root","pdf","1TTLmuRWellResMinQ2_10",false,3,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/LYSO-1TTL-muRwellRes-pythia8Q2-TRonly/output_TRKRS.root","pdf","1TTLmuRWellResMinQ2_10",false,4,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/LYSO-1TTL-muRwellRes-pythia8Q2-TRonly/output_TRKRS.root","pdf","1TTLmuRWellResMinQ2_10",false,5,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/LYSO-1TTL-muRwellRes-pythia8Q2-TRonly/output_TRKRS.root","pdf","1TTLmuRWellResMinQ2_10",false,0,1)'
#     
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/LYSO-1TTL-nomRes-pythia8Q2-TRonly/output_TRKRS.root","pdf","1TTLnomResMinQ2_10",false,0,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/LYSO-1TTL-nomRes-pythia8Q2-TRonly/output_TRKRS.root","pdf","1TTLnomResMinQ2_10",false,3,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/LYSO-1TTL-nomRes-pythia8Q2-TRonly/output_TRKRS.root","pdf","1TTLnomResMinQ2_10",false,4,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/LYSO-1TTL-nomRes-pythia8Q2-TRonly/output_TRKRS.root","pdf","1TTLnomResMinQ2_10",false,5,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/LYSO-1TTL-nomRes-pythia8Q2-TRonly/output_TRKRS.root","pdf","1TTLnomResMinQ2_10",false,0,1)'

#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/2TTL-nomRes-pythia8Q100-TRonly/output_TRKRS.root","pdf","2TTLnomResMinQ2_10",false,0,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/2TTL-nomRes-pythia8Q100-TRonly/output_TRKRS.root","pdf","2TTLnomResMinQ2_10",false,3,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/2TTL-nomRes-pythia8Q100-TRonly/output_TRKRS.root","pdf","2TTLnomResMinQ2_10",false,4,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/2TTL-nomRes-pythia8Q100-TRonly/output_TRKRS.root","pdf","2TTLnomResMinQ2_10",false,5,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/2TTL-nomRes-pythia8Q100-TRonly/output_TRKRS.root","pdf","2TTLnomResMinQ2_10",false,0,1)'

    root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKRS.root","pdf","2TTLnomResMinQ2_1-100",false,0,0)'
    root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKRS.root","pdf","2TTLnomResMinQ2_1-100",false,3,0)'
    root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKRS.root","pdf","2TTLnomResMinQ2_1-100",false,4,0)'
    root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKRS.root","pdf","2TTLnomResMinQ2_1-100",false,5,0)'
    root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKRS.root","pdf","2TTLnomResMinQ2_1-100",false,0,1)'
   
  fi
  if  [ $2 == "effi" ]; then
    basepath=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/LYSO-1TTL-muRwellRes-pythia8Q2-TRonly/output_TRKEFF.root","pdf","1TTLmuRWellResMinQ2_10",0)'
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/LYSO-1TTL-muRwellRes-pythia8Q2-TRonly/output_TRKEFF.root","pdf","1TTLmuRWellResMinQ2_10",1)'
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/LYSO-1TTL-muRwellRes-pythia8Q2-TRonly/output_TRKEFF.root","pdf","1TTLmuRWellResMinQ2_10",2)'
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/LYSO-1TTL-nomRes-pythia8Q2-TRonly/output_TRKEFF.root","pdf","1TTLnomResMinQ2_10",0)'
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/LYSO-1TTL-nomRes-pythia8Q2-TRonly/output_TRKEFF.root","pdf","1TTLnomResMinQ2_10",1)'
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/LYSO-1TTL-nomRes-pythia8Q2-TRonly/output_TRKEFF.root","pdf","1TTLnomResMinQ2_10",2)'
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/2TTL-nomRes-pythia8Q100-TRonly/output_TRKEFF.root","pdf","2TTLnomResMinQ2_10",0)'
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/2TTL-nomRes-pythia8Q100-TRonly/output_TRKEFF.root","pdf","2TTLnomResMinQ2_10",1)'
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/2TTL-nomRes-pythia8Q100-TRonly/output_TRKEFF.root","pdf","2TTLnomResMinQ2_10",2)'
    root -x -l -q -b 'trackingeffi.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKEFF.root","pdf","2TTLnomResMinQ2_1-100",0)'
    root -x -l -q -b 'trackingeffi.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKEFF.root","pdf","2TTLnomResMinQ2_1-100",1)'
    root -x -l -q -b 'trackingeffi.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKEFF.root","pdf","2TTLnomResMinQ2_1-100",2)'
  fi
  
  if [ $2 == "compare" ]; then
#     root -b -x -q -l 'compare_trackingreso_Pythia.C("comparisonTTLLayer.txt","TTLimprovement","pdf",kFALSE)'
#     root -b -x -q -l 'compare_trackingreso_Pythia.C("comparison2TTLLayer.txt","TTLimprovement2TTL","pdf",kFALSE)'
#     root -b -x -q -l 'compare_trackingreso_Pythia.C("comparisonNTTLLayer.txt","TTLNlayers","pdf",kFALSE)'
    root -b -x -q -l 'compare_trackingreso_Pythia.C("comparisonNTTLLayerslim.txt","TTLNlayers2","pdf",kFALSE)'
#     root -b -x -q -l 'compare_trackingeffi.C("comparisonTTLeffNlayer.txt","TTLimprovementNlayers","pdf")'
#     root -b -x -q -l 'compare_trackingeffi.C("comparisonTTLeff.txt","TTLimprovement","pdf")'
#     root -b -x -q -l 'compare_trackingeffi.C("comparison2TTLeff.txt","TTLimprovement2TTL","pdf")'
    
  fi
fi
