#root -x -l -q -b trackingreso.C
#root resolutionJETS.C

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

'
trkrid=2
trksource=0
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-pythia8-10x100-q2-1-to-100-geoOption5/output_TRKRS.root","pdf","pythia-centprod-geo5-q2-1-100", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-pythia8-10x100-q2-1-to-100-geoOption6/output_TRKRS.root","pdf","pythia-centprod-geo6-q2-1-100", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-pythia8-10x100-q2-100-geoOption5/output_TRKRS.root","pdf","pythia-centprod-geo5-q2-100", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-pythia8-10x100-q2-100-geoOption6/output_TRKRS.root","pdf","pythia-centprod-geo6-q2-100", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption5/output_TRKRS.root","pdf","elec-geo5", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-singleElectron-p-0.3-to-20-geoOption6/output_TRKRS.root","pdf","elec-geo6", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption5/output_TRKRS.root","pdf","pion-geo5", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/cades-singlePion-p-0.3-to-20-geoOption6/output_TRKRS.root","pdf","pion-geo6", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/production-singleElectron-p-0-to-20/output_TRKRS.root","pdf","elec-centprod", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/production-singlePion-p-0-to-20/output_TRKRS.root","pdf","pion-centprod", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/production-pythia8-10x100-q2-1-to-100/output_TRKRS.root","pdf","pythia-centprod-q2-1-100", true,'$trkrid','$trksource')'
root -x -l -q -b 'trackingreso_Pythia.C("'$basepath'/production-pythia8-10x100-q2-100/output_TRKRS.root","pdf","pythia-centprod-q2-100", true,'$trkrid','$trksource')'

'

#make effi with other trackset
