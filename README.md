README.md

Project Title
Exploring chronic and transient tumor hypoxia for predicting the efficacy of hypoxia-activated pro-drugs

Description
This code simulates the three-compound treatment of a given schedule. The compound include: the hypoxia-activated pro-drug (HAP) and/or metabolic sensitizer (Sens), and/or vasodilator (Vaso). The user can specify which (if any) of these compounds will be used and the time of their administration. The HAP must be activated in the hypoxic areas within the tumor tissue and only the active form of the drug is lethal. The Sens and Vaso transiently increase the extent of hypoxia within the tissue. As an outcome, the figures showing treatment progression are produced and the number of dead cells is reported. 


The code is generated in Matlab


Getting started
Dependencies
MATLAB software (the R2020b version on a Mac computer was used for all testing)

Installing
Matlab: this code does not need installation

File manifest
The following set of functions are available
i)	runScript.m                -- set up model parameters, calls runHAPSensVaso.m
ii)	runHAPSensVaso.m  -- run treatment schedule, is called by runScript.m 


Executing the program
Matlab: 
Download all MATLAB files (add all files and folders to your MATLAB path)
Run the m-file routines


Authors
Shreya Mathur,
Shannon Chen,
Katarzyna Rejniak; Katarzyna.Rejniak@moffitt.org


Version History
0.1
Initial release


License
This project is licensed under the GNU General Public License v3.0

