inputE=/media/nschmidt/local/AnalysisSoftwareEIC/visualizationGeomAndClusters/plots/ECals/outputPhiEtaECals.root
inputH=/media/nschmidt/local/AnalysisSoftwareEIC/visualizationGeomAndClusters/plots/HCals/outputPhiEtaHCals.root
root -x -l -b -q 'plotGeometry_Paper.C("'$inputE'","'$inputH'")'
