#! /usr/bin/env bash

CHIP=1
CHANNEL="A1"

/au/standa/home
/au/standa/move $(/au/standa/maps/read_map.sh /au/standa/maps/20220206/standa_map_extra.dat $CHIP $CHANNEL)
/au/standa/zero

TAGNAME="reference-$CHANNEL"
FILENAMEX="$TAGNAME.scan_standa_x.dat"
FILENAMEY="$TAGNAME.scan_standa_y.dat"

POS=$(seq -4000 160 4000)

REPEAT=10

rm -rf $FILENAMEX
for X in $POS; do
    /au/standa/move $X 0

    for i in $(seq 1 $REPEAT); do
	OUTPUT=$(/au/measure/pulser_rate.sh $CHIP $CHANNEL)
	echo "position_x = $X $OUTPUT" | tee -a $FILENAMEX
    done
    
done
/au/measure/tree.sh $FILENAMEX position_x counts_on counts_off period_on period_off rate_on rate_off

rm -rf $FILENAMEY
for Y in $POS; do
    /au/standa/move 0 $Y

    for i in $(seq 1 $REPEAT); do
	OUTPUT=$(/au/measure/pulser_rate.sh $CHIP $CHANNEL)
	echo "position_y = $Y $OUTPUT" | tee -a $FILENAMEY
    done
    
done
/au/measure/tree.sh $FILENAMEY position_y counts_on counts_off period_on period_off rate_on rate_off
