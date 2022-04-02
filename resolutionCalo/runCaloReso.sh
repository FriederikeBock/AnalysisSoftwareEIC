clusterizername=MA
basepath=/media/nschmidt/SSD/treeProcessingOutput

if [ $1 = "noteplots" ]; then

    input=SingleParticleNewGeom
    caloname=FEMC-wMat
    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_CRH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'

    caloname=BECAL-wMat
    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_CRH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'

    caloname=EEMC-wMat
    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_CRH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'

    input=HCALStandaloneSingleParticle
    caloname=LFHCAL
    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_CRH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'

    caloname=CHCAL
    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_CRH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'

    caloname=HCALIN
    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_CRH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'

#    input=FOCAL$filled-BH2-FCTB-FCTungsten-SimplePion-NEW
#    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_RH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE,"FCTB-FCTungsten-Pion-NEW")'

#    input=FOCAL$filled-BH2-FCTB-SimplePion-NEW
#    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_RH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE,"FCTB-Pion-NEW")'

fi


if [ $1 = "BEMC_leadglass" ]; then
    caloname=BECAL-wMat-PbGl
    #root -x -l -q -b 'energyResolutionCalorimeters.C("'/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_TTLGEO_7_EEMAPUPDATE_BCLG/output_CRH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'
    caloname=BECAL-wMat-PbGl-NoNCell
    root -x -l -q -b 'energyResolutionCalorimeters.C("'/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_TTLGEO_7_EEMAPUPDATE_BCLG_NONCELL/output_CRH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'
fi


caloname=FEMC
clusterizername=MA
basepath=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing_centralProduction_raymond

if [ $1 = "new" ]; then

    input=production-singlePion-p-0-to-20
    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_RH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'

    input=FOCAL$filled-BH2-FCTB-FCTungsten-SimplePion-NEW
    #root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_RH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE,"FCTB-FCTungsten-Pion-NEW")'

    input=FOCAL$filled-BH2-FCTB-SimplePion-NEW
    #root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_RH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE,"FCTB-Pion-NEW")'

fi

clusterizername=MA
basepath=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing
if [ $1 = "calibration" ]; then
    input=CALIBRATION_ECAL_ELECTRON
    caloname=FEMC
    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_RH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'
    caloname=EEMC
    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_RH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'
    caloname=BECAL
    root -x -l -q -b 'energyResolutionCalorimeters.C("'$basepath/$input/output_RH.root'","'$caloname'","'$clusterizername'","pdf",kTRUE)'
fi
