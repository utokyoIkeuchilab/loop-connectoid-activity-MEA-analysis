# loop-connectoid-activity-MEA-analysis
This repository contains custom codes used for the Multi electrode array (MEA) data analysis for the paper:

Multi-Organoid Loop Cerebral Connectoids Exhibit Enhanced Neuronal Network Dynamics and Sequence-Specific Entrainment

by

Tomoya Duenki (1,2,3,4) , Yoshiho Ikeuchi (1,2,3,4,*)

1 Institute of Industrial Science, The University of Tokyo, Meguro, Tokyo 153-8505, Japan. 2 Institute for AI and Beyond, The University of Tokyo, Bunkyo, Tokyo 113-8655, Japan. 3 Department of Chemistry and Biotechnology, The University of Tokyo, Bunkyo, Tokyo 113-8655, Japan. 4 LIMMS, CNRS-Institute of Industrial Science, IRL 2820, The University of Tokyo, Tokyo, Japan.


It contains the following 3 main scripts:
- Raw trace analysis (Main_raw_trace_analysis.m):
    - Preprocessing of data
    - Spike detection
    - Plotting of Signals
    - Time frequency analysis
    - Wavelet coherence
    - Frequency separation
- Correlation analysis (Main_correlation_analysis.m):
    - Summerizing of spiking data within organoid
    - Inter organoid correlation (direct & indirect)
    - Intra organoid correlation
    - Plotting spiking activity and correlation matrix
- Burst analysis (Main_burst_analysis.m):
    - Burst detection
    - Burst frequency
    - Interburst interval
    - Interburst interval coefficient of variation
    - Burst peak size
    - Burst size coefficient of variation
    - Consecutive bursts
    - Spontaneous activity transient (SAT) Duration

  
Example data set:

Example data can be found in the "Example data" folder which contains a .mat file with 10 min raw signal traces from 4 organoids that were connected and differentiated for 12 weeks ('example_signal_trace.mat'). It also contains a .mat file with extracted spiking activity from 64 electrodes of the same recording (example_spike_data.mat).


All analysis were performed using MATLAB 2020a.

Required Toolboxes:
-  Signal Processing Toolbox
-  Wavelet Toolbox
-  Statistics and Machine Learning Toolbox
-  Parallel Computing Toolbox
-  Image Processing Toolbox