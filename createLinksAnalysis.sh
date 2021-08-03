EICANA=""

if [ $1 = "fbock" ]; then
  EICANA="/home/fbock/eic/developEIC/AnalysisSoftwareEIC"
fi

mkdir trackingresolution/
mkdir treeAnalysis/
mkdir resolutionHCAL/
mkdir clusterProperties/
mkdir common/
mkdir plotEtaPhi/
ln -s $EICANA/trackingresolution/*.C trackingresolution/
# ln -s $EICANA/trackingresolution/*.h trackingresolution/
ln -s $EICANA/trackingresolution/*.sh trackingresolution/
ln -s $EICANA/treeAnalysis/*.h treeAnalysis/
ln -s $EICANA/treeAnalysis/*.cxx treeAnalysis/
ln -s $EICANA/treeAnalysis/*.C treeAnalysis/
ln -s $EICANA/treeAnalysis/*.sh treeAnalysis/
ln -s $EICANA/clusterProperties/*.C clusterProperties/
ln -s $EICANA/plotEtaPhi/*.C plotEtaPhi/
ln -s $EICANA/plotEtaPhi/*.sh plotEtaPhi/
ln -s $EICANA/resolutionHCAL/*.C resolutionHCAL/
ln -s $EICANA/common/*.h common/
ln -s $EICANA/*.sh .
