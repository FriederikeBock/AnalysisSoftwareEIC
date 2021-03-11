EICANA=""

if [ $1 = "fbock" ]; then
  EICANA="/home/fbock/eic/developEIC/AnalysisSoftwareEIC"
fi

mkdir trackingresolution/
mkdir treeAnalysis/
mkdir common/
ln -s $EICANA/trackingresolution/*.C trackingresolution/
# ln -s $EICANA/trackingresolution/*.h trackingresolution/
ln -s $EICANA/trackingresolution/*.sh trackingresolution/
ln -s $EICANA/treeAnalysis/*.h treeAnalysis/
ln -s $EICANA/treeAnalysis/*.cxx treeAnalysis/
ln -s $EICANA/treeAnalysis/*.C treeAnalysis/
ln -s $EICANA/treeAnalysis/*.sh treeAnalysis/
ln -s $EICANA/common/*.h common/
