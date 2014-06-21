## SJS. shell script to vary beta (boltzmann) for submission a qsub script.

SCRIPT=siminf.qsub

# establish initial condition variable
START=1
STOP=6


# submit initial condition
qsub $SCRIPT

# submit rest
for ((i=$START;i<$STOP;i++)); do
    NEXT=$[i+1]
    sed -i "s/BETA=$i/BETA=$NEXT/g" $SCRIPT
    qsub $SCRIPT
done

# Reset
sed -i "s/BETA=$NEXT/BETA=$START/g" $SCRIPT


