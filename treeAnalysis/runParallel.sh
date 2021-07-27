#! /bin/bash
# NUMTHREADS=16
# INPUT="/home/tristan/ecce/Singularity/analysis/AnalysisSoftwareEIC/treeAnalysis/test_prod.txt"

# invoke with ./runParalle.sh -j NUMBER_OF_THREADS -f LIST_OF_FILES

echo "IMPORTANT NOTE: Only tested for jet resolution studies, other studies may not combined histograms in a sensible way."
echo "If verified for other studies, please leave a note here."

# Stop all child processes when an appropriate signal is received
function stop {
    for PID in "${PIDS[@]}"; do
        echo Killing $PID
        kill $PID
    done
    exit
}
trap stop SIGHUP SIGINT SIGTERM

while getopts j:f: option; do   # Get core count and file list
    case "${option}" in
        j) NUMTHREADS=${OPTARG};;
        f) INPUT=$(readlink -f ${OPTARG});;     # Exapand to absolute path so changing directories later doesn't break things
    esac
done

if [ -z ${NUMTHREADS+x} ] || [ -z ${INPUT+x} ]; then
    echo "Invoke with $0 -j NUMBER_OF_THREADS -f LIST_OF_FILES"
    exit
fi

# Staging areas
mkdir -p treeProcessing/input_files
mkdir -p treeProcessing/logs

# Split the input into NUMTHREADS sections
cd treeProcessing/input_files
split -d --number=l/$NUMTHREADS $INPUT
cd ../..

PIDS=()

echo ""
echo "Starting ${NUMTHREADS} processes..."
echo "Progress can be monitored with 'tail -f treeProcessing/logs/*.log'"
# Start NUMTHREADS instances of the analyzer 
for (( i=0; i<$NUMTHREADS; i++ )); do
  printf -v J "%02d" $i
  root -b -x -q -l 'treeProcessing.C("treeProcessing/input_files/x'$J'", "geometry.root", "'$J'", true, true)' > "treeProcessing/logs/$J.log" 2> "treeProcessing/logs/$J.err" &
  PIDS+=( $! )
done;

wait    # Wait for all root processes to return
echo "Done processing, merging histograms..."

# Merge histograms with hadd
for FILE in  output_CS.root output_HITS.root output_JRH.root output_TMSTUD.root output_TRKRS.root output_CLSIZER.root output_EventObservables.root output_JetObservables.root output_RH.root output_TRKEFF.root output_TRKS_Comparison.root; do
    FILEPATHS=()
    for (( i=0; i<$NUMTHREADS; i++)); do 
        printf -v FILEPATH "treeProcessing/%02d/%s" $i $FILE
        FILEPATHS+=( $FILEPATH )
    done
    echo ${FILEPATHS[@]}
    hadd -j $NUMTHREADS -f treeProcessing/$FILE ${FILEPATHS[@]}
done

echo ""
echo "Done"

# Ryzen 3700X, 4,000,000 events from SIDIS high Q^2
# Cores | Time mm:ss
#     1 |
#     2 |
#     4 |
#     8 |
#    16 | 100:12
