input=/media/nschmidt/local/AnalysisSoftwareEIC/clusterstudies/plotsEoverP/2022_04_04/output_EoverP_ExampleBins.root
root -x -l -b -q 'eoverpplotting_paper.C("'$input'","pdf")'
