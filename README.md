

## Format of the input files

- for log files:
write them: log-HV1-HV2-...-HVn.lvm
in a folder logFiles/$detectorName/date

- for MCA files:
write them: spectrum-HV1-HV2-...-HVn-CoarseGain.mca
in a folder MCA/$detectorName/date


# Use of this analysis code

- Check detector information in _DetectorInfo.C (in particular, geometry, calibration constant...)

- If Calibration needs to be done: ``` ./GetCalibration.sh ```

- Then type ```./Process.sh```and follow the instructions
