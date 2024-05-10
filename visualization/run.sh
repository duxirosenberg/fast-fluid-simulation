#run script simply 
NX=8
NY=8
NZ=8
DIRECTION=27
STREAM=4
COLLISON=1
LBM=1
MOMENTUM=1
RESET_DATAFILE=0



# compile in build folder ...

# ./timing $DIRECTION $NX $NY $NZ $STREAM $MOMENTUM $COLLISON $LBM $RESET_DATAFILE
# ./build/timing $DIRECTION $NX $NY $NZ $STREAM $MOMENTUM 
# for NXYZ in 2 4 8 16; do ./timing $DIRECTION $NXYZ $NXYZ $NXYZ $STREAM $MOMENTUM $COLLISON $LBM $RESET_DATAFILE; done

# ./timing 9 1 1 1 4 1 1 1 1
# reset datafile ..
../build/timing 9 1 1 1 0 0 0 0 1
for NXYZ in 2 4 8 16 32 64;
do
    for DIR in 9 15 27
    do 
        echo "currently doing $DIR and $NXYZ"
        ../build/timing $DIR $NXYZ $NXYZ $NXYZ $STREAM $MOMENTUM $COLLISON $LBM $RESET_DATAFILE
    done
done
