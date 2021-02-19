#!/bin/bash



#detector="ZZBOT"
#detector="RD3SP3"
#detector="RD3SP1"
#detector="LittleChinese"
#detector="RD3SP5"
detector="MGEM1"

# Start here

root -l -q "CalibrationOrtec.C(\"$detector\")"

