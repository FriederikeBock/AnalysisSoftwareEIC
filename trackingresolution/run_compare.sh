#root -x -l -q -b trackingreso.C
#root resolutionJETS.C
trkrid=2
#effi
#root -x -l -q -b 'compare_trackingeffi.C("/media/nschmidt/local/AnalysisSoftwareEIC/trackingresolution/comparePion.txt","pions","pdf")'

#root -x -l -q -b 'compare_trackingeffi.C("/media/nschmidt/local/AnalysisSoftwareEIC/trackingresolution/compareElectron.txt","electrons","pdf")'
root -x -l -q -b 'compare_trackingeffi.C("/media/nschmidt/local/AnalysisSoftwareEIC/trackingresolution/comparePythia.txt","pythia","pdf")'


#reso
trksrc=0
#root -x -l -q -b 'compare_trackingreso_Pythia.C("/media/nschmidt/local/AnalysisSoftwareEIC/trackingresolution/comparePionReso.txt","pions","pdf",true)'


#root -x -l -q -b 'compare_trackingreso_Pythia.C("/media/nschmidt/local/AnalysisSoftwareEIC/trackingresolution/compareElectronReso.txt","electrons","pdf",true)'
#root -x -l -q -b 'compare_trackingreso_Pythia.C("/media/nschmidt/local/AnalysisSoftwareEIC/trackingresolution/comparePythiaReso.txt","pythia","pdf",true)'
