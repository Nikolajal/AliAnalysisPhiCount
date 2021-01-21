#! /bin/bash

./bash/MCG_PhiQuickReference.sh
mkdir result_MCG_PhiQuickReference || exit 1
mkdir result_MCG_PhiQuickReference/logs

strun=0
nruns=1
njobs=5
nevents=100000

echo "[INFO] Starting production"
        

for run in $(seq $strun $(($strun + $nruns - 1))); do
    
    for option in {-1,0,1,2,3,4,5,6,7}; do

    ### wait if there are too many jobs running
    while true; do
        bkgjobs=$(jobs | grep MCG_PhiQuickReference | wc -l | xargs)
        if [ $bkgjobs -lt $njobs ]; then
            break
        fi
        echo "[INFO] sleep while waiting for a free job slot"
        sleep 60
    done
    
    runid=$(printf "%09d" $(( run + 1000000000 * (option +1) )))
    seed=$((123456789 + $run * 2))
    
    echo "[INFO] Starting run $runid with production option $option"

    touch result_MCG_PhiQuickReference/logs/MCG_PhiQuickReference_$runid.log
    ./exe/MCG_PhiQuickReference result_MCG_PhiQuickReference/MCG_PhiQuickReference_$runid $nevents $seed $option >& result_MCG_PhiQuickReference/logs/MCG_PhiQuickReference_$runid.log &

    sleep 1s

    done
done

exit 0
