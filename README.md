# Explore ARISE controller injection amounts to climate modes

## Files
**controller/commonrutines.py** : Common routines used in the controller 

**controller/driver.py** : Controller driver

**controller/main_controller.py** : Main controller document. Calls other controller scripts (commonrutines.py, driver.py, PIcontrol.py). Allows a user to give the controller a map of temperature (K) and returns injection amounts at the 4 locations.

**controller/PIcontrol.py** : Controller calculations can be found in this script. Creates the .txt output file.

**createControllerInputs/CalculateBasestate.py** : Calculates and saves the Basestate (Smoothed climate chnage) data.

**createControllerInputs/CalcualteENSO.py** : Calculated ENSO from temperature anomalies

**createContollerInputs/main_input.py** : Main file used to add temperature anomalies to basestates and then feed temperature map to main_controller.py.

