# MDD participants compared to healthy age- and gender-matched controls on the n-back working memory task

INQUISIT:
Inquisit software was used to run the experiment. 
Tasks were presented in the following order:
  2-back task practice;
  Eyes-closed resting state;
  IAPS;
  Eyes-open resting state;
  0-, 1-, and 2-back tasks;
  

EEG PROCESSING PIPELINE:
Scripts can be run using MATLAB. 
Although not tested, open-source sofware Octave may also be able to run the scripts.

The functioning of each of the scripts is as follows:

MDDvsCTRL_pipeline: calls all other scripts/functions and coordinates processing steps.
MDDvsCTRL_readTMSi: reads EEG data from the TMSi Polybench (.Poly5) format using the 'tms_read' function.
MDDvsCTRL_preprocess: preprocesses the data by removing line noise (at 50Hz), and implements a bandpass filter. 
MDDvsCTRL_events: identifies trigger events in the data using a marker channel, e.g. onsent of resting-state
MDDvsCTRL_trials: segments the continuous EEG data into epochs. For task-based EEG (IAPS and n-back task), the epochs are equal to the length of the trials. For resting-state EEG, the epochs are one second in length.
MDDvsCTRL_trialrejection: semi-automated rejection of trials based on data range and z-scores
MDDvsCTRL_ICA...: four scripts are provided (ICAdata, ICAauto, ICAplot, ICAcheck) which, respectively, concatenate epoched data for ICA analysis, suggest identified components for rejection, plot components for visual inspection, and finally coordinate the action of the three previously described scripts. 


SPSS STATISTICAL ANALYSES:
Mixed-effects repeated measures model and regression analyses were performed using data in a longform data format. 


QUERIES:
stevan.nikolin@unsw.edu.au
stevan.nikolin@gmail.com
