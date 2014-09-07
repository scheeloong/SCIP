#!/bin/sh

# Initialize variables to be randomly generated 
# Note: IN bash, need n=20 NOT n = 20, whitespace matters
YRANGE=20
ZRANGE=20
NRANGE=25
XRANGE=30
n=$RANDOM
Xlb=0
Xub=0
let "n %=$NRANGE"
Ylb=$RANDOM
let "Ylb %=$YRANGE" 
Yub=$RANDOM
let "Yub %=$YRANGE" 
let "Yub +=$Ylb"

Zlb=$RANDOM
let "Zlb %=$ZRANGE" 
Zub=$RANDOM
let "Zub %=$ZRANGE" 
let "Zub +=$Zlb"

# Initialize to an empty file 
echo "$n $Ylb $Yub $Zlb $Zub" > ElementTestRandomGenerated.ns
count=1
while [ "$count" -le $n ]      # Generate n random integers bounds
do
    Xlb=$RANDOM
    let "Xlb %=$XRANGE" 
    Yub=$RANDOM
    let "Xub %=$XRANGE" 
    let "Xub +=$Xlb"
    echo "$Xlb $Xub" >> ElementTestRandomGenerated.ns
    let "count += 1"
done
