# -*- coding: utf-8 -*-
"""

Main script to read and analyse data for correlation analysis propulsion-gait speed

Sint Maartenskliniek study ID: 1028_OpenGo

Last update: September 2023
Authors: C.J. Ensink, J. Biere (v.09-2023), S. van Outvorst (v. 07-2023), C.J. Ensink (initial version)


"""

# Import dependencies
import pandas as pd
import glob
import os
import numpy as np
from readmarkerdata import readmarkerdata
from gaiteventdetection import gaiteventdetection
from gaitcharacteristics import spatiotemporals, propulsion

# Create initial dictionaries to save the marker data, forceplate data, gait event detection, spatiotemporals
markerdata = dict()
analogdata = dict()
gaitevents = dict()
gaitcharacteristics = dict()

# # # # # HEALTHY PARTICIPANTS # # # # #

# Location of the healthy participants information file
healthy_participant_information_path = 'V:/research_stage/2023SamvanOutvorst/VII_Matlab/MOCAP_analysis/Subject information gezond.xlsx'

# Location of the optical motion capture (.c3d) files
healthy_data_path = 'V:/research_stage/2023SamvanOutvorst/VII_Matlab/Python Propulsie/c3d_gezond'
# General file name structure
trial_name_structure = "pp*_S*.c3d"

# Import particpant information
healthy_participant_information = pd.read_excel(healthy_participant_information_path)
healthy_participants = healthy_participant_information["ID"].tolist()

# Create initial dictionaries to save the average propulsion and step length asymmetry per particpant, per trial
propulsion_average_impulse = dict()
propulsion_average_peak = dict()
braking_average_impulse = dict()
braking_average_peak = dict()
step_length_average_asymmetry = dict()

# Analyze the data
for file in glob.glob(os.path.join(healthy_data_path, trial_name_structure)):
    # Get participant ID and trial name from the .c3d filename
    filename_without_extension = os.path.splitext(os.path.basename(file))[0]
    participantID, trialname = filename_without_extension.split("_")
    
    healthyfilename = 'healthy_'+filename_without_extension
    
    # Print
    print(f"Analyze data of file: {file}")

    # Get participant bodyweight from healthy participant information file
    bodyweight = healthy_participant_information[healthy_participant_information["ID"] == participantID]["Weight"].values[0]
    
    # Read markerdata, read force plate data
    markerdata[filename_without_extension], fs_markerdata, analogdata[filename_without_extension], fs_analogdata = readmarkerdata(file, analogdata=True)
    # Identify gait events
    gaitevents[filename_without_extension] = gaiteventdetection(markerdata[filename_without_extension], fs_markerdata, algorithmtype='velocity', trialtype='treadmill', debugplot=True)
    # Calculate spatiotemporal gait parameters
    gaitcharacteristics[filename_without_extension] = spatiotemporals(markerdata[filename_without_extension], gaitevents[filename_without_extension], debugplot=True)
    # Calculate propulsion
    gaitevents[filename_without_extension], gaitcharacteristics[filename_without_extension], analogdata[filename_without_extension] = propulsion(gaitevents[filename_without_extension], gaitcharacteristics[filename_without_extension], analogdata[filename_without_extension], bodyweight=bodyweight, debugplot=True, plot_title=healthyfilename)

    # Average propolsion (impulse) of the left and right side of each trial
    propulsion_impulse_left_avg = np.nanmean(gaitcharacteristics[filename_without_extension]['Propulsion left'][:, 2])
    propulsion_impulse_right_avg = np.nanmean(gaitcharacteristics[filename_without_extension]['Propulsion right'][:, 2])
    try:
        propulsion_average_impulse[participantID][trialname] = np.mean([propulsion_impulse_left_avg, propulsion_impulse_right_avg])
    except KeyError:
        try:
            propulsion_average_impulse[participantID] = dict()
            propulsion_average_impulse[participantID][trialname] = np.mean([propulsion_impulse_left_avg, propulsion_impulse_right_avg])
        except:
            print('Something went wrong in saving the average propulsion impulse')
    
    # Average propulsion (peak) of the let and right side of each trial
    propulsion_peak_left_avg = np.nanmean(gaitcharacteristics[filename_without_extension]['Peak propulsion left'][:,1])
    propulsion_peak_right_avg = np.nanmean(gaitcharacteristics[filename_without_extension]['Peak propulsion right'][:,1])
    try:
        propulsion_average_peak[participantID][trialname] = np.mean([propulsion_peak_left_avg, propulsion_peak_right_avg])
    except KeyError:
        try:
            propulsion_average_peak[participantID] = dict()
            propulsion_average_peak[participantID][trialname] = np.mean([propulsion_peak_left_avg, propulsion_peak_right_avg])
        except:
            print('Something went wrong in saving the average propulsion peak')
    
    # Average braking (impulse) of the left and right side of each trial
    braking_left_avg = np.nanmean(gaitcharacteristics[filename_without_extension]['Braking left'][:, 2])
    braking_right_avg = np.nanmean(gaitcharacteristics[filename_without_extension]['Braking right'][:, 2])
    try:
        braking_average_impulse[participantID][trialname] = np.mean([braking_left_avg, braking_right_avg])
    except KeyError:
        try:
            braking_average_impulse[participantID] = dict()
            braking_average_impulse[participantID][trialname] = np.mean([braking_left_avg, braking_right_avg])
        except:
            print('Something went wrong in saving the average braking impulse')
    
    # Average braking (peak) of the let and right side of each trial
    braking_peak_left_avg = np.nanmean(gaitcharacteristics[filename_without_extension]['Peak braking left'][:,1])
    braking_peak_right_avg = np.nanmean(gaitcharacteristics[filename_without_extension]['Peak braking right'][:,1])
    try:
        braking_average_peak[participantID][trialname] = np.mean([braking_peak_left_avg, braking_peak_right_avg])
    except KeyError:
        try:
            braking_average_peak[participantID] = dict()
            braking_average_peak[participantID][trialname] = np.mean([braking_peak_left_avg, braking_peak_right_avg])
        except:
            print('Something went wrong in saving the average peak braking')

# Convert dictionary to pandas DataFrame
df_propulsion_average_impulse = pd.DataFrame(dict((k, propulsion_average_impulse[k]) for k in healthy_participants if k in propulsion_average_impulse))
df_braking_average_impulse = pd.DataFrame(dict((k, braking_average_impulse[k]) for k in healthy_participants if k in braking_average_impulse))
df_propulsion_average_peak = pd.DataFrame(dict((k, propulsion_average_peak[k]) for k in healthy_participants if k in propulsion_average_peak))
df_braking_average_peak = pd.DataFrame(dict((k, braking_average_peak[k]) for k in healthy_participants if k in braking_average_peak))

# Set excel output filenames
filename1 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Healthy_propulsion_impulse.xlsx'
filename2 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Healthy_braking_impulse.xlsx'
filename3 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Healthy_propulsion_peak.xlsx'
filename4 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Healthy_braking_peak.xlsx'

# Export DataFrames to Excel
df_propulsion_average_impulse.to_excel(filename1, index=True)
df_braking_average_impulse.to_excel(filename2, index=True)
df_propulsion_average_peak.to_excel(filename3, index=True)
df_braking_average_peak.to_excel(filename4, index=True)





# # # # # STROKE PARTICIPANTS # # # # #

# Location of the stroke participants information file
stroke_participant_information_path = 'V:/research_stage/2023SamvanOutvorst/VII_Matlab/MOCAP_analysis/Subject information.xlsx'

# Location of the optical motion capture (.c3d) files
stroke_data_path = 'V:/research_stage/2023SamvanOutvorst/VII_Matlab/MOCAP_analysis/VICON_data'
# General file name structure
trial_name_structure = "pp*_S*.c3d"

# Import particpant information
stroke_participant_information = pd.read_excel(stroke_participant_information_path)
stroke_participants = list()

# Analyze the data
for file in glob.glob(os.path.join(stroke_data_path, trial_name_structure)):
    # Get participant ID and trial name from the .c3d filename
    filename_without_extension = os.path.splitext(os.path.basename(file))[0]
    participantID, trialname = filename_without_extension.split("_")
    
    strokefilename = 'stroke_'+filename_without_extension
    strokeID = 'stroke_'+ participantID
    stroke_participants.append(strokeID)
    
    # Print
    print(f"Analyze data of file: {file}")

    # Get participant bodyweight from stroke participant information file
    bodyweight = stroke_participant_information[stroke_participant_information["ID"] == participantID]["Weight"].values[0]
    
    # Read markerdata, read force plate data
    markerdata[strokefilename], fs_markerdata, analogdata[strokefilename], fs_analogdata = readmarkerdata(file, analogdata=True)
    # Identify gait events
    gaitevents[strokefilename] = gaiteventdetection(markerdata[strokefilename], fs_markerdata, algorithmtype='velocity', trialtype='treadmill', debugplot=True)
    # Calculate spatiotemporal gait parameters
    gaitcharacteristics[strokefilename] = spatiotemporals(markerdata[strokefilename], gaitevents[strokefilename])
    # Calculate propulsion
    try:
        gaitevents[strokefilename], gaitcharacteristics[strokefilename], analogdata[strokefilename] = propulsion(gaitevents[strokefilename], gaitcharacteristics[strokefilename], analogdata[strokefilename], bodyweight=bodyweight, debugplot=True, plot_title=strokefilename)
    except:
        print(f'Error in trial:  {strokefilename}')
    
    # Spatiotemporeel
    # Steplength 
    Steplength_L_avg = np.nanmean(gaitcharacteristics[strokefilename]['Steplength left (mm)'][:,1])
    Steplength_R_avg = np.nanmean(gaitcharacteristics[strokefilename]['Steplength right (mm)'][:,1])
    try:
        Steplength_avg[strokeID][trialname] = dict()
        Steplength_avg[strokeID][trialname]['left'] = Steplength_L_avg
        Steplength_avg[strokeID][trialname]['right'] = Steplength_R_avg
    except KeyError:
        try:
            Steplength_avg[strokeID] = dict()
            Steplength_avg[strokeID][trialname] = dict()
            Steplength_avg[strokeID][trialname]['left'] = Steplength_L_avg
            Steplength_avg[strokeID][trialname]['right'] = Steplength_R_avg       
        except:
            print('Something went wrong in saving the average steplength')   
            
        
        
    Steplength_left[strokeID][trialname] = dict()
    Steplength_left[strokeID][trialname]['Steplength left (mm)'] = (gaitcharacteristics['stroke_pp09_S11']['Steplength left (mm)'][:,1])
    
    
    
    (gaitcharacteristics['stroke_pp09_S11']['Steplength left (mm)'][:,1])
    
    
    # Average propolsion (impulse) of the left and right side of each trial
    propulsion_impulse_left_avg = np.nanmean(gaitcharacteristics[strokefilename]['Propulsion left'][:, 2])
    propulsion_impulse_right_avg = np.nanmean(gaitcharacteristics[strokefilename]['Propulsion right'][:, 2])
    try:
        propulsion_average_impulse[strokeID][trialname] = dict()
        propulsion_average_impulse[strokeID][trialname]['left'] = propulsion_impulse_left_avg
        propulsion_average_impulse[strokeID][trialname]['right'] = propulsion_impulse_right_avg
    except KeyError:
        try:
            propulsion_average_impulse[strokeID] = dict()
            propulsion_average_impulse[strokeID][trialname] = dict()
            propulsion_average_impulse[strokeID][trialname]['left'] = propulsion_impulse_left_avg
            propulsion_average_impulse[strokeID][trialname]['right'] = propulsion_impulse_right_avg
        except:
            print('Something went wrong in saving the average propulsion impulse')
    
    # Average propulsion (peak) of the left and right side of each trial
    propulsion_peak_left_avg = np.nanmean(gaitcharacteristics[strokefilename]['Peak propulsion left'][:,1])
    propulsion_peak_right_avg = np.nanmean(gaitcharacteristics[strokefilename]['Peak propulsion right'][:,1])
    try:
        propulsion_average_peak[strokeID][trialname] = dict()
        propulsion_average_peak[strokeID][trialname]['left'] = propulsion_peak_left_avg
        propulsion_average_peak[strokeID][trialname]['right'] = propulsion_peak_right_avg
    except KeyError:
        try:
            propulsion_average_peak[strokeID] = dict()
            propulsion_average_peak[strokeID][trialname] = dict()
            propulsion_average_peak[strokeID][trialname]['left'] = propulsion_peak_left_avg
            propulsion_average_peak[strokeID][trialname]['right'] = propulsion_peak_right_avg
        except:
            print('Something went wrong in saving the average propulsion peak')
    
    # Average braking (impulse) of the left and right side of each trial
    braking_left_avg = np.nanmean(gaitcharacteristics[strokefilename]['Braking left'][:, 2])
    braking_right_avg = np.nanmean(gaitcharacteristics[strokefilename]['Braking right'][:, 2])
    try:
        braking_average_impulse[strokeID][trialname] = dict()
        braking_average_impulse[strokeID][trialname]['left'] = braking_left_avg
        braking_average_impulse[strokeID][trialname]['right'] = braking_right_avg
    except KeyError:
        try:
            braking_average_impulse[strokeID] = dict()
            braking_average_impulse[strokeID][trialname] = dict()
            braking_average_impulse[strokeID][trialname]['left'] = braking_left_avg
            braking_average_impulse[strokeID][trialname]['right'] = braking_right_avg
        except:
            print('Something went wrong in saving the average braking impulse')
    # Average braking (peak) of the let and right side of each trial
    braking_peak_left_avg = np.nanmean(gaitcharacteristics[strokefilename]['Peak braking left'][:,1])
    braking_peak_right_avg = np.nanmean(gaitcharacteristics[strokefilename]['Peak braking right'][:,1])
    try:
        braking_average_peak[strokeID][trialname] = dict()
        braking_average_peak[strokeID][trialname]['left'] = braking_peak_left_avg
        braking_average_peak[strokeID][trialname]['right'] = braking_peak_right_avg
    except KeyError:
        try:
            braking_average_peak[strokeID] = dict()
            braking_average_peak[strokeID][trialname] = dict()
            braking_average_peak[strokeID][trialname]['left'] = braking_peak_left_avg
            braking_average_peak[strokeID][trialname]['right'] = braking_peak_right_avg
        except:
            print('Something went wrong in saving the average braking')

# Convert dictionary to pandas DataFrame for left side outcome measures
df_dict = dict()
for p in stroke_participants:
    df_dict[p] = dict()
    for trial in propulsion_average_impulse[p]:
        df_dict[p][trial] = propulsion_average_impulse[p][trial]['left']
df_propulsion_average_impulse = pd.DataFrame(df_dict)
df_dict = dict()
for p in stroke_participants:
    df_dict[p] = dict()
    for trial in braking_average_impulse[p]:
        df_dict[p][trial] = braking_average_impulse[p][trial]['left']
df_braking_average_impulse = pd.DataFrame(df_dict)
df_dict = dict()
for p in stroke_participants:
    df_dict[p] = dict()
    for trial in propulsion_average_peak[p]:
        df_dict[p][trial] = propulsion_average_peak[p][trial]['left']
df_propulsion_average_peak = pd.DataFrame(df_dict)
df_dict = dict()
for p in stroke_participants:
    df_dict[p] = dict()
    for trial in braking_average_peak[p]:
        df_dict[p][trial] = braking_average_peak[p][trial]['left']
df_braking_average_peak = pd.DataFrame(df_dict)

# Set excel output filenames for left outcome measures
filename1 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Stroke_propulsion_impulse_left.xlsx'
filename2 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Stroke_braking_impulse_left.xlsx'
filename3 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Stroke_propulsion_peak_left.xlsx'
filename4 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Stroke_braking_peak_left.xlsx'

# Export DataFrames to Excel for left outcome measures
df_propulsion_average_impulse.to_excel(filename1, index=True)
df_braking_average_impulse.to_excel(filename2, index=True)
df_propulsion_average_peak.to_excel(filename3, index=True)
df_braking_average_peak.to_excel(filename4, index=True)

# Convert dictionary to pandas DataFrame for right side outcome measures
df_dict = dict()
for p in stroke_participants:
    df_dict[p] = dict()
    for trial in propulsion_average_impulse[p]:
        df_dict[p][trial] = propulsion_average_impulse[p][trial]['right']
df_propulsion_average_impulse = pd.DataFrame(df_dict)
df_dict = dict()
for p in stroke_participants:
    df_dict[p] = dict()
    for trial in braking_average_impulse[p]:
        df_dict[p][trial] = braking_average_impulse[p][trial]['right']
df_braking_average_impulse = pd.DataFrame(df_dict)
df_dict = dict()
for p in stroke_participants:
    df_dict[p] = dict()
    for trial in propulsion_average_peak[p]:
        df_dict[p][trial] = propulsion_average_peak[p][trial]['right']
df_propulsion_average_peak = pd.DataFrame(df_dict)
df_dict = dict()
for p in stroke_participants:
    df_dict[p] = dict()
    for trial in braking_average_peak[p]:
        df_dict[p][trial] = braking_average_peak[p][trial]['right']
df_braking_average_peak = pd.DataFrame(df_dict)

# Set excel output filenames for right outcome measures
filename5 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Stroke_propulsion_impulse_right.xlsx'
filename6 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Stroke_braking_impulse_right.xlsx'
filename7 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Stroke_propulsion_peak_right.xlsx'
filename8 = 'V:/research_reva_studies/1028_OpenGo/II_Onderzoeksdata/Python/Stroke_braking_peak_right.xlsx'

# Export DataFrames to Excel for right outcome measures
df_propulsion_average_impulse.to_excel(filename5, index=True)
df_braking_average_impulse.to_excel(filename6, index=True)
df_propulsion_average_peak.to_excel(filename7, index=True)
df_braking_average_peak.to_excel(filename8, index=True)


