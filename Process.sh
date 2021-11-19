#!/bin/bash

clean()
{
filesToDelete="*.so *.d *.pcm *ACLiC* *.in *.out"
for entry in $filesToDelete
do
if test -f $entry;
then
#echo Deleting $entry
rm $entry
fi
done
}


testMca() {
for file in *mca; do
	#echo $file
	if test -f $file;
	then found=1
	fi
done
}

diplayGoodDir() {
mainDir=$PWD
cd MCA/$detector
mcaDir=$PWD
for folder in `find ./  -type d  | sort -f`; do
folder="${folder:3}"
cd $folder
found=0
testMca
if [ $found = 1 ]; then
echo $folder
fi
cd $mcaDir
done
cd $mainDir
}

detectorList="LittleChinese MGEM1 MGEM3 RD3SP1 RD3SP3 ZZBOT"

# Start here

detectorOK=0
if ! [ -z $1 ]
then
for entry in $detectorList
do
if [ $1 = $entry ]
then
detectorOK=1
fi
done
fi

if [ $detectorOK != 1 ]
then
echo "Please type ./Process.sh 'detector name' 'folder' where detector can be:"
for entry in $detectorList
do
echo $entry
done
exit
fi

detector=$1

echo
echo "Detector: $detector"
echo

if [ -z $2 ]
then
echo "Please type ./Process.sh 'detector name' 'folder' where folder can be:"
diplayGoodDir
exit
fi

clean
folder=$2

figureFolder="Figures/$detector/$folder"
if [ ! -d $figureFolder ];
then
mkdir $figureFolder
fi

if [ ! -d $figureFolder ];
then
echo "Could not create the folder to store the figures"
exit
fi


root -l -q "Analyse.C(\"$detector\", \"$folder\")"
clean


