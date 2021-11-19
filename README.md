

## Format of the input files

- for log files currents files, deprecated unless done with Labview)):
write them: log-HV1-HV2-...-HVn.lvm
in a folder logFiles/$detectorName/date

- for current files with the python script:
pico-HV1-HV2-...-with.csv and pico-HV1-HV2-...-without.csv (with and without source)

- for MCA files:
write them: spectrum-HV1-HV2-...-HVn-CoarseGain.mca
in a folder MCA/$detectorName/date


# Using this code

- Check detector information in Include/DetectorInfo.C (in particular, geometry, calibration constant...)

- If Calibration needs to be done: ``` ./GetCalibration.sh ```

- Then type ```./Process.sh``` and follow the instructions
