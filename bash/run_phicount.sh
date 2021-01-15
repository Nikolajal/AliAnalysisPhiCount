./GeneratorMC.sh
mkdir MCgen || exit 1

strun=0
nruns=100
njobs=12
nevents=1000000

for run in $(seq $strun $(($strun + $nruns - 1))); do

    ### wait if there are too many jobs running
    while true; do
        bkgjobs=$(jobs | grep MCG_PhiAnalysis | wc -l | xargs)
        if [ $bkgjobs -lt $njobs ]; then
            break
        fi
        echo "[---] sleep while waiting for a free job slot"
        sleep 60
    done
    
    runid=$(printf "%05d" $run)
    seed=$((123456789 + $run * 2))
    
    echo "[---] starting run: $runid"

    ./exe/MCG_PhiAnalysis result/outGeneratorMC_$runid.root $nevents $runid >& /dev/null &

    sleep 1s

done
exit 0
