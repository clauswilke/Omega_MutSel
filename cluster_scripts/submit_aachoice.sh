## SJS. shell script to vary numaa from 4->12. starts with NUMAA=4

# establish initial condition variable
start=4

# submit initial NUMAA=4 condition
qsub aachoice.qsub

# submit lots
for ((i=4;i<12;i++)); do
    next=$[i+1]
    sed -i "s/NUMAA=$i/NUMAA=$next/g" aachoice.qsub 
    qsub aachoice.qsub
done

# Reset
sed -i "s/NUMAA=$next/NUMAA=$start/g" aachoice.qsub


