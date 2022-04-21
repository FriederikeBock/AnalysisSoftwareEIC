if [ $1 == "run" ] ; then
  suffix="pdf";
datep=2022_04_08
  #inputfile=/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/ECal_inputs.txt
    echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionEEMC-wMatFitting/$datep/ResolutionPostProcessing.root EEMC EEMC 0\n/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionBECAL-wMatFitting/$datep/ResolutionPostProcessing.root BECAL BEMC 1\n/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionFEMC-wMatFitting/$datep/ResolutionPostProcessing.root FEMC FEMC 2" > ECal_inputs.txt
  root -b -x -l -q 'energyResolutionCalorimeters_comparePlot.C("ECal_inputs.txt","ECals","'$suffix'")'
#fi

  #inputfile=/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/HCal_inputs.txt
    echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionCHCALFitting/$datep/ResolutionPostProcessing.root HCALOUT OHCAL 1\n/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionLFHCALFitting/$datep/ResolutionPostProcessing.root LFHCAL LFHCAL 2" > HCal_inputs.txt
  root -b -x -l -q 'energyResolutionCalorimeters_comparePlot.C("HCal_inputs.txt","HCals","'$suffix'")'


  #inputfile=/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/ECal_inputs2.txt
    echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionBECAL-wMat-PbGlFitting/ResolutionPostProcessing.root BECAL BEMC_PbGl 0\n/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionBECAL-wMatFitting/ResolutionPostProcessing.root BECAL BEMC_SciGl 1" > ECal_inputs2.txt
  #root -b -x -l -q 'energyResolutionCalorimeters_comparePlot.C("ECal_inputs2.txt","BEMC_comp","'$suffix'")'

elif [ $1 == "comp" ] ; then
    #echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionBECAL-wMat-PbGlFitting/ResolutionPostProcessing.root BEMC_PbGl\n/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionBECAL-wMatFitting/ResolutionPostProcessing.root BEMC_SciGl" > ECal_inputs_comp.txt
    echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionBECAL-wMat-PbGl-NoNCellFitting/ResolutionPostProcessing.root BEMC_PbGl\n/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionBECAL-wMatFitting/ResolutionPostProcessing.root BEMC_SciGl" > ECal_inputs_comp.txt

    root -b -x -l -q 'compare_EnergyResol.C("ECal_inputs_comp.txt","PbGlCompNoNCell","BECAL","png",3)'
    root -b -x -l -q 'compare_EnergyResol.C("ECal_inputs_comp.txt","PbGlCompNoNCell","BECAL","png",14)'
    root -b -x -l -q 'compare_EnergyResol.C("ECal_inputs_comp.txt","PbGlCompNoNCell","BECAL","png",23)'


elif [ $1 == "papercomp" ] ; then
    #echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionBECAL-wMat-PbGlFitting/ResolutionPostProcessing.root BEMC_PbGl\n/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionBECAL-wMatFitting/ResolutionPostProcessing.root BEMC_SciGl" > ECal_inputs_comp.txt
    echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionEEMC-wMatFitting/ResolutionPostProcessing.root EEMC 10" > ECal_inputs_comp.txt
    echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionBECAL-wMatFitting/ResolutionPostProcessing.root BEMC 8" >> ECal_inputs_comp.txt
    echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionFEMC-wMatFitting/ResolutionPostProcessing.root FEMC 20" >> ECal_inputs_comp.txt
    echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionCHCALFitting/ResolutionPostProcessing.root OHCAL 8" >> ECal_inputs_comp.txt
    echo -e "/media/nschmidt/local/AnalysisSoftwareEIC/resolutionCalo/plotsResolutionLFHCALFitting/ResolutionPostProcessing.root LFHCAL 22" >> ECal_inputs_comp.txt

   root -b -x -l -q 'compare_EnergyResol_paper.C("ECal_inputs_comp.txt","PaperPlots","","pdf")'
   # root -b -x -l -q 'compare_EnergyResol.C("ECal_inputs_comp.txt","PbGlCompNoNCell","BECAL","png",14)'
   # root -b -x -l -q 'compare_EnergyResol.C("ECal_inputs_comp.txt","PbGlCompNoNCell","BECAL","png",23)'

fi
