#!/bin/bash

##----------------------------------------------------------------------------##
## Authors: Nick Brazeau
## This script is used get rid of barcodes
##----------------------------------------------------------------------------##

if [ $# -eq 0 ]; then
  echo "Usage: symlinker [OPTIONS]"
  echo "Try symlinker --help"
  exit 1
fi

##----------------------------------------------------------------------------##
## Collect inputs
##----------------------------------------------------------------------------#

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -h|--help)
    HELP=true
    ;;
    -I|--input)
    INPUT="$2"
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

##----------------------------------------------------------------------------##
## Help file
##----------------------------------------------------------------------------##


if [ "$HELP" = true ]; then

	echo "Usage: symlinker [OPTIONS]"
	echo ""
	echo "  -I, --input		A dir for splitting on underscores (_)."
	exit 0

fi



##----------------------------------------------------------------------------##
## write out symlinks
##----------------------------------------------------------------------------##
len=`ls $INPUT`
for i in $len;
  do
    end=`echo $i | awk -F'_L001_' '{print $2}'`
    smpl=`echo $i | awk -F'[A|C|T|G]' '{print $1}'`
    smpl=`echo $smpl | sed 's/_/-/g'`
    out=`echo ${smpl}${end} | sed 's/-R/_R/'`
    echo "renaming sample from=$INPUT/$i to=$INPUT/$out"
    mv $INPUT/$i $INPUT/$out

done
