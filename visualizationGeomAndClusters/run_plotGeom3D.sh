input=/media/nschmidt/SSD/simulationOutputCADES/tst/output_TTLGEO_7_HITS_LFTAILC_BCNG_EEMAPUPDATE_fixed_SimplePiZero/geometry.root
detector=ECals
caloindex=-1
root -b -x -l -q 'plotGeometry3D.C("'$input'","pdf","'$detector'",'$caloindex')'

detector=HCals
caloindex=-2
root -b -x -l -q 'plotGeometry3D.C("'$input'","pdf","'$detector'",'$caloindex')'

