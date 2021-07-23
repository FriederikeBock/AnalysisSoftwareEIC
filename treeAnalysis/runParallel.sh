#! /bin/bash
NUMTHREADS=8
INPUT="/home/tristan/ecce/Singularity/analysis/AnalysisSoftwareEIC/treeAnalysis/test_prod.txt"

# Stop all child processes when an appropriate signal is received
function stop {
    for PID in "${PIDS[@]}"; do
        echo Killing $PID
        kill $PID
        exit
    done
}
trap stop SIGHUP SIGINT SIGTERM

# Staging areas
mkdir -p treeProcessing/input_files
mkdir -p treeProcessing/logs

# Split the input into NUMTHREADS sections
cd treeProcessing/input_files
split -d --number=l/$NUMTHREADS $INPUT
cd ../..

PIDS=()

# Start NUMTHREADS instances of the analyzer 
for (( i = 0; i < $NUMTHREADS; i++ )); do
  printf -v J "%02d" $i
  root -b -x -q -l 'treeProcessing.C("treeProcessing/input_files/x'$J'", "geometry.root", "'$J'", true, true)' > "treeProcessing/logs/$J.log" 2> "treeProcessing/logs/$J.err" &
  PIDS+=( $! )
done;

wait    # Wait for all root processes to return
echo "Done processing, merging histograms..."

# Merge histograms with hadd
for FILE in  output_CS.root output_HITS.root output_JRH.root output_TMSTUD.root output_TRKRS.root output_CLSIZER.root output_EventObservables.root output_JetObservables.root output_RH.root output_TRKEFF.root output_TRKS_Comparison.root; do
    FILEPATHS=()
    for (( i = 0; i < $NUMTHREADS; i++)); do 
        printf -v FILEPATH "treeProcessing/%02d/%s" $i $FILE
        FILEPATHS+=( $FILEPATH )
    done
    # echo ${FILEPATHS[@]}
    hadd -j $NUMTHREADS -f treeProcessing/$FILE ${FILEPATHS[@]}
done

echo "Done"


# Cores | Time mm:ss
#     1 | 42:01
#     2 |
#     4 |
#     7 | 16:46
#     8 |
