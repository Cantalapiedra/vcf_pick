inputfile=$1
parent1=$2
parent2=$3
allowhet=$4

BASEDIR=$(dirname $0)

$BASEDIR/src/check_polymorph.py $1 $2 $3 $4 > $1".polymorph."$2"_"$3"_"$4.out 2> $1".polymorph."$2"_"$3"_"$4.err

## END
