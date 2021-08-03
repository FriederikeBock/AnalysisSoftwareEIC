# for i in {0..99};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/G4EICDetector_eventtree.root","pdf","SinglePion",10,'$i')'
# done

# for i in {0..99};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_SinglePion_30GeV/G4EICDetector_eventtree.root","pdf","SinglePion",30,'$i')'
# done

# for i in {0..99};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_SinglePion_2GeV/G4EICDetector_eventtree.root","pdf","SinglePion",2,'$i')'
# done
# 
for i in {0..99};
do 
  root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_SingleProton_2GeV/G4EICDetector_eventtree.root","pdf","SingleProton",2,'$i')'
done

for i in {0..99};
do 
  root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_SingleProton_30GeV/G4EICDetector_eventtree.root","pdf","SingleProton",30,'$i')'
done

# for i in {0..99};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_SingleProton_100GeV/G4EICDetector_eventtree.root","pdf","SingleProton",100,'$i')'
# done
