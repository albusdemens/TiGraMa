# Script to make separate folders for the images relative to different projections

for (( i=0; i <= 48; i++ ))
do
  mkdir "/Isolated_blobs_$i"
done

# Separate the files in different folders
cd /data/alcer/Data_analysis/Reconstruction_Murofi/Fe_2015/Isolated_blobs_100/

for (( i=0; i <= 9; i++ ))
do
mv "Isolated_blob_00$i"* "/Isolated_blobs_$i"
done
for (( i=10; i <= 48; i++ ))
do
mv "Isolated_blob_0$i"* "/Isolated_blobs_$i"
done
