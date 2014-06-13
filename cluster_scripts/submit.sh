## SJS. shell script to vary numaa for submission a qsub script.

SCRIPT=convergence.qsub

# establish initial condition variable
START=4
STOP=15


# submit initial condition
qsub $SCRIPT

# submit rest
for ((i=$START;i<$STOP;i++)); do
    NEXT=$[i+1]
    sed -i "s/NUMAA=$i/NUMAA=$NEXT/g" $SCRIPT
    qsub $SCRIPT
done

# Reset
sed -i "s/NUMAA=$NEXT/NUMAA=$START/g" $SCRIPT


