#inputFileTR="/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_SINGLEMULTIPION_TTLGEO7_Match/output_TRKEFF.root"
#inputFileCS="/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_SINGLEMULTIPION_TTLGEO7_Match/output_CS.root"
inputFileTR="/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_SINGLEMULTIPIONELEC_TTLGEO7_Match/output_TRKEFF.root"
inputFileCS="/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_SINGLEMULTIPIONELEC_TTLGEO7_Match/output_CS.root"

  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","LFHCAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALOUT","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALIN","pdf","SinglePart","",0,0)'
#  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EHCAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","FEMC","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","BECAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EEMC","pdf","SinglePart","",0,0)'
