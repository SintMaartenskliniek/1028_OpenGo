Code used for data-analysis of the TZO project, part of the ZonMw-TopZorg (TZO) grant nr. 10070022010004 

** FULL ARTICLE UNDER SUBMISSION: Joost Biere, Brenda E. Groen, Carmen J. Ensink, Jorik Nonnekes, Noël L.W. Keijsers: Gait speed-dependent modulation of paretic versus non-paretic propulsion in persons with chronic stroke

Main script to read and analyse optical motion capture data (VICON based)
Possibility to include analog (forceplate) data

Main script to read and analyse the anteroposterior vector of the ground reaction force; propulsion peak and propulsion pulse during gait

INPUT 
Motion capture (VICON based) and analog (forceplate) data of healthy and/or stroke participants 
Posibility for multiple conditions, specify in trial name structure 
.c3d file format required

OUTPUT
Excel files with propulsion pulse and propulsion peak per leg, per condition, per participant

**********************************************************************************************************************************************

Run the main_propulsion.py code to analyse the data 
	- The main scripts requires the eventdetection.py, gaitcharacteristics.py and readmarkerdata.py scripts to run
	- The folder exampleDATA contains two folders with example data of a healthy and a stroke subject, walking on a number of gait speeds
	- The folder exampleDATA contains two excel files 'Subject information healthy' and 'Subject information stroke', needed in de main script

Explanation of the code is commented through the scripts. Explanation of the example data is provided in the ReadMe file in the folder "exampleDATA".