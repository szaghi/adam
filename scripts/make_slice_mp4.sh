#!/bin/bash

if [ $# -eq 0 ] ; then
  echo "Usage"
  echo "      ./make_slice_mp4.sh xy"
  echo "options:"
  echo "      xy is the (zero padded) number of slice"
  echo ""
  echo "example:"
  echo "      ./make_slice_mp4.sh 01"
  exit
fi

sn=$1
sd='slice_'$sn'_png/'
mkdir -p $sd
for f in *slice_$sn*dat
do
   echo "make slice png of $f"
   bf=`basename $f .dat`
   cat slice_png_template.gp | sed -e "s/BNAME/$bf/g" > $bf.gp
   gnuplot $bf.gp
   mv $bf.png $sd
   rm -f $bf.gp
done
cd $sd
ffmpeg -framerate 10 -pattern_type glob -i "*.png" slice_$sn.mp4
mv slice_$sn.mp4 ../
cd ..
rm -rf $sd
