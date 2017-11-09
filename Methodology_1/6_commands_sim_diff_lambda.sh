# Script to make separate folders for the images relative to different projections

for (( i=0; i <= 48; i++ ))
do
  mkdir "Comb_blobs_ordered/Comb_blobs_$i"
done

# Separate the files in different folders
cd Comb_blobs/

for (( i=0; i <= 9; i++ ))
do
cp "Comb_blob_00$i"* "Comb_blobs_ordered/Comb_blobs_$i/"
done
for (( i=10; i <= 48; i++ ))
do
cp "Comb_blob_0$i"* "Comb_blobs_ordered/Comb_blobs_$i/"
done
