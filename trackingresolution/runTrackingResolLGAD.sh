
if [ $1 = "plots" ]; then
  fileLBLLGAD=rootfiles/histosTrackingResol_LBLwithLGAD.root
  fileLBL=rootfiles/histosTrackingResol_defaultLBL.root
  fileFST=rootfiles/histosTrackingResol_defaultFST.root
  fileFSTLGAD=rootfiles/histosTrackingResol_FST-LGAD.root
        
  root -b -x -q -l trackingreso_flatPt.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Hist\"\,kFALSE\)
  root -b -x -q -l trackingreso_flatPt.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Fit\"\)

  root -b -x -q -l trackingreso_flatPt.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist\"\,kFALSE\)
  root -b -x -q -l trackingreso_flatPt.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Fit\"\)

  root -b -x -q -l trackingreso_flatPt.C\(\"$fileFST\"\,\"pdf\"\,\"defaultFST-Hist\"\,kFALSE\)
  root -b -x -q -l trackingreso_flatPt.C\(\"$fileFST\"\,\"pdf\"\,\"defaultFST-Fit\"\)

  root -b -x -q -l trackingreso_flatPt.C\(\"$fileFSTLGAD\"\,\"pdf\"\,\"FSTwithLGAD-Hist\"\,kFALSE\)
  root -b -x -q -l trackingreso_flatPt.C\(\"$fileFSTLGAD\"\,\"pdf\"\,\"FSTwithLGAD-Fit\"\)
  
elif [ $1 = "tree" ]; then
  
  root -x -b -q -l analyseTreeForTrackingResol.C\(\"defaultLBL\"\,\"\"\,\"filesLBL.txt\"\)
  root -x -b -q -l analyseTreeForTrackingResol.C\(\"defaultFST\"\,\"\"\,\"filesFST.txt\"\)
  root -x -b -q -l analyseTreeForTrackingResol.C\(\"LBLwithLGAD\"\,\"\"\,\"filesLGAD.txt\"\)
  root -x -b -q -l analyseTreeForTrackingResol.C\(\"FST-LGAD\"\,\"\"\,\"filesFST-FTTL.txt\"\)

fi
