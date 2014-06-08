## SJS. shell script to vary numaa from 2 -> 10. starts with NUMAA=2

qsub aachoice.qsub
sed -i "s/NUMAA=2/NUMAA=3/g" aachoice.qsub
qsub aachoice.qsub
sed -i "s/NUMAA=3/NUMAA=4/g" aachoice.qsub
qsub aachoice.qsub
sed -i "s/NUMAA=4/NUMAA=5/g" aachoice.qsub
qsub aachoice.qsub
sed -i "s/NUMAA=5/NUMAA=6/g" aachoice.qsub
qsub aachoice.qsub
sed -i "s/NUMAA=6/NUMAA=7/g" aachoice.qsub
qsub aachoice.qsub
sed -i "s/NUMAA=7/NUMAA=8/g" aachoice.qsub
qsub aachoice.qsub
sed -i "s/NUMAA=8/NUMAA=9/g" aachoice.qsub
qsub aachoice.qsub
sed -i "s/NUMAA=9/NUMAA=10/g" aachoice.qsub
qsub aachoice.qsub

# Reset
sed -i "s/NUMAA=10/NUMAA=2/g" aachoice.qsub


