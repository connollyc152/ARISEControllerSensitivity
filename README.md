# Quantifying the impact of internal variability on the CESM control algorithm for stratospheric aerosol injection

The original code for the Assessing Responses and Impacts of Solar (ARISE-SAI) controller can be found at: https://doi.org/10.5281/zenodo.6471092

## Files
**controller/commonrutines.py** : Common routines used in the controller 

**controller/driver.py** : Controller driver

**controller/main_controller.py** : Main controller document. Calls other controller scripts (commonrutines.py, driver.py, PIcontrol.py). Allows a user to give the controller a map of temperature (K) and returns injection amounts at the 4 locations.

**controller/PIcontrol.py** : Controller calculations can be found in this script. Creates the .txt output file.

**createControllerInputs/CalculateBasestate.py** : Calculates and saves the Basestate (Smoothed climate chnage) data.

**createControllerInputs/CalcualteENSO.py** : Calculate temperature anomalies associated with the ENSO

**createControllerInputs/CalcualteSAM.py** : Calculate temperature anomalies associated with the SAM

**createControllerInputs/CalcualteNAO.py** : Calculate temperature anomalies associated with the NAO

**createControllerInputs/CalcualteVolc.py** : Calculated temperature anomalies associated with a the Pinatubo eruption

**createContollerInputs/main_input.py** : Main file used to add temperature anomalies to basestates and then feed temperature map to main_controller.py.

