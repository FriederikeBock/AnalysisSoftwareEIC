#root -x -l -q -b trackingreso.C
#root resolutionJETS.C
if [ $1 == "Nico" ]; then
basepath=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing_centralProduction_raymond
trkrid=0

#root -x -l -q -b 'trackingeffi.C("'$basepath'/production-singleElectron-p-0-to-20/output_TRKEFF.root","pdf","elec-centprod",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/production-singlePion-p-0-to-20/output_TRKEFF.root","pdf","pion-centprod",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption5/output_TRKEFF.root","pdf","elec-geo5",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption6/output_TRKEFF.root","pdf","elec-geo6",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption5/output_TRKEFF.root","pdf","pion-geo5",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption6/output_TRKEFF.root","pdf","pion-geo6",'$trkrid')'

root -x -l -q -b 'trackingeffi.C("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_SINGLEMULTIPIONELEC_TTLGEO7_Match/output_TRKEFF.root","pdf","singlepionelec_projfix",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_SINGLEMULTIPION_TTLGEO7_Match/output_TRKEFF.root","pdf","singleelec_projfix",'$trkrid')'

#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-pythia8-10x100-q2-1-to-100-geoOption5/output_TRKEFF.root","pdf","pythia-centprod-geo5-q2-1-100",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-pythia8-10x100-q2-1-to-100-geoOption6/output_TRKEFF.root","pdf","pythia-centprod-geo6-q2-1-100",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-pythia8-10x100-q2-100-geoOption5/output_TRKEFF.root","pdf","pythia-centprod-geo5-q2-100",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-pythia8-10x100-q2-100-geoOption6/output_TRKEFF.root","pdf","pythia-centprod-geo6-q2-100",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption5/output_TRKEFF.root","pdf","elec-geo5",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption6/output_TRKEFF.root","pdf","elec-geo6",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption5/output_TRKEFF.root","pdf","pion-geo5",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption6/output_TRKEFF.root","pdf","pion-geo6",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/production-singleElectron-p-0-to-20/output_TRKEFF.root","pdf","elec-centprod",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/production-singlePion-p-0-to-20/output_TRKEFF.root","pdf","pion-centprod",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/production-pythia8-10x100-q2-1-to-100/output_TRKEFF.root","pdf","pythia-centprod-q2-1-100",'$trkrid')'
#root -x -l -q -b 'trackingeffi.C("'$basepath'/production-pythia8-10x100-q2-100/output_TRKEFF.root","pdf","pythia-centprod-q2-100",'$trkrid')'

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

if [ $1 == "Fredi" ]; then
  if [ $2 == "reso" ]; then
    basepath=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing

#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/SingleElectronEEMCPCarbonEOP/output_TRKRS.root","pdf","SingleElectron",false,0,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/SinglePionEEMCPCarbonEOP/output_TRKRS.root","pdf","SinglePion",false,0,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/SingleElectronEEMCPCarbonEOP/output_TRKRS.root","pdf","SingleElectron",true,0,0)'
#     root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/SinglePionEEMCPCarbonEOP/output_TRKRS.root","pdf","SinglePion",true,0,0)'

    root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/SingleParticle_EOP_PbGlasG4/output_TRKRS.root","pdf","SingleElectron",true,0,0)'
    root -b -x -q -l 'trackingreso_Pythia.C("'$basepath'/SingleParticle_EOP_PbGlasG4/output_TRKRS.root","pdf","SinglePion",true,0,0)'
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
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKEFF.root","pdf","2TTLnomResMinQ2_1-100",0)'
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKEFF.root","pdf","2TTLnomResMinQ2_1-100",1)'
#     root -x -l -q -b 'trackingeffi.C("'$basepath'/2TTL-nomRes-pythia8Q1-100-TRonly/output_TRKEFF.root","pdf","2TTLnomResMinQ2_1-100",2)'

    root -x -l -q -b 'trackingeffi.C("'$basepath'/SingleParticle_EOP_PbGlasG4/output_TRKEFF.root","pdf","1TTLMB_1-100",0)'

  fi
  
  if [ $2 == "compare" ]; then
#     root -b -x -q -l 'compare_trackingreso_Pythia.C("comparisonTTLLayer.txt","TTLimprovement","pdf",kFALSE)'
#     root -b -x -q -l 'compare_trackingreso_Pythia.C("comparison2TTLLayer.txt","TTLimprovement2TTL","pdf",kFALSE)'
#     root -b -x -q -l 'compare_trackingreso_Pythia.C("comparisonNTTLLayer.txt","TTLNlayers","pdf",kFALSE)'
#     root -b -x -q -l 'compare_trackingreso_Pythia.C("comparisonNTTLLayerslim.txt","TTLNlayers2","pdf",kFALSE)'

    #fit
    root -b -x -q -l 'compare_trackingreso_Pythia.C("comparison2TTLLayer.txt","TTLimprovement2TTL","pdf",kTRUE)'
    root -b -x -q -l 'compare_trackingreso_Pythia.C("comparisonNTTLLayer.txt","TTLNlayers","pdf",kTRUE)'
    root -b -x -q -l 'compare_trackingreso_Pythia.C("comparisonNTTLLayerslim.txt","TTLNlayers2","pdf",kTRUE)'

    
#     root -b -x -q -l 'compare_trackingeffi.C("comparisonTTLeffNlayer.txt","TTLimprovementNlayers","pdf")'
#     root -b -x -q -l 'compare_trackingeffi.C("comparisonTTLeff.txt","TTLimprovement","pdf")'
#     root -b -x -q -l 'compare_trackingeffi.C("comparison2TTLeff.txt","TTLimprovement2TTL","pdf")'
    
  fi
fi
