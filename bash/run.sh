#! /bin/bash

./bash/GeneratorMC.sh
mkdir -p result/MC_Production/

strun=971
nruns=$((2500-$strun))
njobs=3
nevents=100000

for run in $(seq $strun $(($strun + $nruns - 1))); do

    ### wait if there are too many jobs running
    while true; do
        bkgjobs=$(ps aux | grep MCG_PhiAnalysis | wc -l | xargs)
        if [ $bkgjobs -lt $(($njobs +1)) ]; then
            break
        fi
        echo "[---] sleep while waiting for a free job slot"
        sleep 60
    done
    
    runid=$(printf "%05d" $run)
    seed=$((123456789 + $run * 2))
    
    echo "[---] starting run: $runid"

    ./exe/MCG_PhiAnalysis result/MC_Production/outGeneratorMC_$runid.root $nevents $runid 0 >& ./result/MC_Production/log.$runid.log &

    sleep 1s

done
exit 0
