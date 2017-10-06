DLDIR=/data/projects/mziemann/geo2mx/firefox_download
WD=/data/projects/mziemann/geo2mx_project/v1/firefox_update

OLD=$WD/old
FILE=sra_result.csv

#export DISPLAY=:0
#xhost localhost

mv $WD/ecoli_$FILE $OLD
rm $DLDIR/$FILE
firefox &
firefox "imacros://run/?m=ecoli.iim" &
sleep 60
mv $DLDIR/$FILE $WD/ecoli_$FILE

exit


rm $DLDIR/$FILE
firefox "imacros://run/?m=athaliana.iim" &
sleep 60
mv $DLDIR/$FILE $DLDIR/athaliana_$FILE

rm $DLDIR/$FILE
firefox "imacros://run/?m=celegans.iim" &
sleep 60
mv $DLDIR/$FILE $DLDIR/celegans_$FILE

rm $DLDIR/$FILE
firefox "imacros://run/?m=dmelanogaster.iim" &
sleep 60
mv $DLDIR/$FILE $DLDIR/dmelanogaster_$FILE

rm $DLDIR/$FILE
firefox "imacros://run/?m=drerio.iim" &
sleep 60
mv $DLDIR/$FILE $DLDIR/drerio_$FILE

rm $DLDIR/$FILE
firefox "imacros://run/?m=ecoli.iim" &
sleep 60
mv $DLDIR/$FILE $DLDIR/ecoli_$FILE

rm $DLDIR/$FILE
firefox "imacros://run/?m=hsapiens.iim" &
sleep 120
mv $DLDIR/$FILE $DLDIR/hsapiens_$FILE

rm $DLDIR/$FILE
firefox "imacros://run/?m=mmusculus.iim" &
sleep 120
mv $DLDIR/$FILE $DLDIR/mmusculus_$FILE

rm $DLDIR/$FILE
firefox "imacros://run/?m=rnovegicus.iim" &
sleep 60
mv $DLDIR/$FILE $DLDIR/rnovegicus_$FILE

rm $DLDIR/$FILE
firefox "imacros://run/?m=scerevisiae.iim" &
sleep 20
mv $DLDIR/$FILE $DLDIR/scerevisiae_$FILE

exit
