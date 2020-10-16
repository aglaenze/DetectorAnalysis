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
for folder in `find ./  -type d`; do
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


#detector="ZZBOT"
#detector="RD3SP3"
#detector="RD3SP1"
#detector="LittleChinese"
detector="RD3SP5"

# Start here

echo
echo "Detector: $detector"
echo

if [ -z $1 ]
then
echo "Please type ./Process.sh $folderName where folder can be:"
diplayGoodDir
exit
fi

clean
folder=$1

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


