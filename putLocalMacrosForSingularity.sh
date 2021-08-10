EICMACROS=""
SPHENIXMACROS=""
SINGULARIYMACROS=""
SINGULARIYCALIB=""
EICCALIB=""

if [ $1 = "fbock" ]; then
  EICMACROS="/home/fbock/eic/developEIC/EICmacros/common"
  SPHENIXMACROS="/home/fbock/eic/developEIC/sPHENIXmacros/common"
  SINGULARIYMACROS="/home/fbock/eic/Singularity/cvmfs/eic.opensciencegrid.org/gcc-8.3/release/release_new/new/rootmacros"
  SINGULARIYCALIB="/home/fbock/eic/Singularity/cvmfs/eic.opensciencegrid.org/gcc-8.3/release/new/share/calibrations"
  EICCALIB="/home/fbock/eic/developEIC/calibrations_eic"
fi

rm $SINGULARIYMACROS/*.C
rm $SINGULARIYCALIB/Beam/*.gdml 
if [ $2 = "link" ]; then
  ln -s $EICMACROS/*.C $SINGULARIYMACROS/
  ln -s $SPHENIXMACROS/*.C $SINGULARIYMACROS/
  ln -s $EICCALIB/Beam/*.gdml $SINGULARIYCALIB/Beam/
elif [ $2 = "copy" ]; then
  cp $EICMACROS/*.C $SINGULARIYMACROS/
  cp $SPHENIXMACROS/*.C $SINGULARIYMACROS/
fi
