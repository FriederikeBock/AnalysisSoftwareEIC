EICMACROS=""
SPHENIXMACROS=""
SINGULARIYMACROS=""

if [ $1 = "fbock" ]; then
  EICMACROS="/home/fbock/eic/developEIC/EICmacros/common"
  SPHENIXMACROS="/home/fbock/eic/developEIC/sPHENIXmacros/common"
  SINGULARIYMACROS="/home/fbock/eic/Singularity/cvmfs/eic.opensciencegrid.org/gcc-8.3/release/release_new/new/rootmacros"
fi

rm $SINGULARIYMACROS/*.C
if [ $2 = "link" ]; then
  ln -s $EICMACROS/*.C $SINGULARIYMACROS/
  ln -s $SPHENIXMACROS/*.C $SINGULARIYMACROS/
elif [ $2 = "copy" ]; then
  cp $EICMACROS/*.C $SINGULARIYMACROS/
  cp $SPHENIXMACROS/*.C $SINGULARIYMACROS/
fi
