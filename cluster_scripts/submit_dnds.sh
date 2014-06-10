## SJS. shell script to vary numaa and submit dnds.qsub

# establish initial condition variable
START=3
STOP=20


# submit initial condition
qsub dnds.qsub

# submit rest
for ((i=$START;i<$STOP;i++)); do
    NEXT=$[i+1]
    sed -i "s/NUMAA=$i/NUMAA=$NEXT/g" dnds.qsub 
    qsub dnds.qsub
done

# Reset
sed -i "s/NUMAA=$NEXT/NUMAA=$START/g" dnds.qsub


