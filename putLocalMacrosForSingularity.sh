EICMACROS=""
SPHENIXMACROS=""
SINGULARIYMACROS=""
SINGULARIYCALIB=""
EICCALIB=""

if [ $1 = "fbockECCE" ]; then
  EICMACROS="/home/fbock/EIC/Software/eccemacros/common"
  SPHENIXMACROS="/home/fbock/EIC/Software/sPHENIXmacros/common"
  SINGULARIYMACROS="/home/fbock/EIC/Singularity_ECCE/cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/release/release_new/new/rootmacros"
  SINGULARIYCALIB="/home/fbock/EIC/Singularity_ECCE/cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/release/new/share/calibrations"
  EICCALIB="/home/fbock/EIC/Software/eccecalibrations"
elif [ $1 = "fbockEIC" ]; then
  EICMACROS="/home/fbock/EIC/Software/eicmacros/common"
  SPHENIXMACROS="/home/fbock/EIC/Software/sPHENIXmacros/common"
  SINGULARIYMACROS="/home/fbock/EIC/Singularity_EIC/cvmfs/eic.opensciencegrid.org/gcc-8.3/release/release_new/new/rootmacros"
  SINGULARIYCALIB="/home/fbock/EIC/Singularity_EIC/cvmfs/eic.opensciencegrid.org/gcc-8.3/release/new/share/calibrations"
  EICCALIB="/home/fbock/EIC/Software/eiccalibrations"
fi

ls $SINGULARIYMACROS/*.C
rm $SINGULARIYMACROS/*.C
# rm $SINGULARIYCALIB/Beam/*.gdml 
rm $SINGULARIYCALIB/BarrelEcal/mapping/*.txt
if [ $2 = "link" ]; then
  ln -s $EICMACROS/*.C $SINGULARIYMACROS/
  if [[ "$1" == *"EIC"* ]]; then
    ln -s $SPHENIXMACROS/*.C $SINGULARIYMACROS/
  fi
#   ln -s $EICCALIB/Beam/*.gdml $SINGULARIYCALIB/Beam/
  ln -s $EICCALIB/BarrelEcal/mapping/*.txt $SINGULARIYCALIB/BarrelEcal/mapping/
elif [ $2 = "copy" ]; then
  cp $EICMACROS/*.C $SINGULARIYMACROS/
  if [[ "$1" == *"EIC"* ]]; then
    cp $SPHENIXMACROS/*.C $SINGULARIYMACROS/
  fi
fi
