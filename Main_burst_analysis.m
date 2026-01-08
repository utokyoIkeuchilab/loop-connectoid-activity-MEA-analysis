%% ============================================================
%  Main burst analysis script
%  ============================================================

clear all; close all; clc;



%% ============================================================
%  Set parameters and import data
%  ============================================================

% Base parameters
total_duration  = 1200; % Total recording duration (seconds)
bin_win         = 100; % Binning window (ms) for spike counts
smoothingfactor = 3;  % Gaussian smoothing strength
bin_number      = total_duration*1000/bin_win; % Number of bins for the full recording
time            = [total_duration/bin_number:total_duration/bin_number:total_duration];

% Parameters for burst calculations
burstthresholdpercent    = 0.1; % Initial relative threshold
continued_IBI_threshold  = 15;  % Max inter-burst interval to be classified as continued bursts (seconds)

% Define electrode location
TL = [1,2,3,4,9,10,11,12,17,18,19,20,25,26,27,28]; % Electrodes of top left organoids
TR = [5,6,7,8,13,14,15,16,21,22,23,24,29,30,31,32]; % Electrodes of top right organoids
BL = [33,34,35,36,41,42,43,44,49,50,51,52,57,58,59,60]; % Electrodes of bottom left organoids
BR = [37,38,39,40,45,46,47,48,53,54,55,56,61,62,63,64]; % Electrodes of bottom right organoids


% Load MATLAB .mat file containing the spike timing data for all 64 electrodes (can be downloaded on Github)
All_spikes = struct2array(load('C:\Users\Tomoya\Projects\papers\loop-connectoid-activity-MEA-analysis\example_data\example_spike_data.mat'));




%% ============================================================
%  Select organoid and pool spikes
%  ============================================================

% Select which organoid to be analyzed
% TL for top left, TR for top right, BL for bottom left and BR for bottom right
selorg = [TL]; 

Spikes = [];
for i = 1:16 % write all spike timings of selected organoid into one vector
    Spikes = [Spikes; All_spikes{selorg(i), 1}];
end

% Convert spike times into binned spike counts
[N,~] = histcounts(Spikes, bin_number);
% Smooth population activity using Gaussian kernel
N2 = conv(N,gausswin(smoothingfactor*(1/(bin_win/1000))));
% Define absolute burst threshold as fraction of maximum activity
burstthreshold = burstthresholdpercent*max(N2);

% Detect population bursts as peaks above threshold
[burst_spikes, burst_locs] = findpeaks(N2,'MinPeakHeight', burstthreshold);



%% ============================================================
%  Compute different parameters based on detected burst location
%  ============================================================

% Mean burst frequency (bursts per second)
burstparameters.meanburstfrequency = size(burst_locs,2)/total_duration;
burstparameters.burstlocation = rot90(burst_locs/(1/(bin_win/1000)),3);

% Inter-burst intervals (seconds)
burstparameters.meanibi = mean(diff(burst_locs/(1/(bin_win/1000))));
burstparameters.ibi = rot90(diff(burst_locs/(1/(bin_win/1000))),3);
burstparameters.ibi_CV = nanstd([burstparameters.ibi; burst_locs(1)/10; total_duration-burst_locs(end)/10])/nanmean([burstparameters.ibi; burst_locs(1)/10;total_duration-burst_locs(end)/10]);

% Burst size
burstparameters.burstspikenumber = rot90(burst_spikes,3);
burstparameters.meanburstsize = mean(burst_spikes);
burstparameters.burstsize_cv = nanstd(burstparameters.burstspikenumber)/nanmean(burstparameters.burstspikenumber);

% Continued burst 
isConsec = burstparameters.ibi  <= continued_IBI_threshold; % Identify consecutive bursts (IBI ≤ threshold)
d          = diff([0; isConsec; 0]); % Find start/end indices of consecutive runs 
runStarts  = find(d == 1); 
runEnds    = find(d == -1) - 1;
burstparameters.consecutiveBurst_lengths = runEnds - runStarts + 2; % Number of bursts per run
burstparameters.consecutiveBurst_durations = burstparameters.burstlocation(runEnds + 1) - burstparameters.burstlocation(runStarts); % Duration of each consecutive burst episode (s)
burstparameters.mean_consecutiveBurst_lengths = mean(burstparameters.consecutiveBurst_lengths); % Mean number of bursts per run
burstparameters.mean_consecutiveBurst_durations = mean(burstparameters.consecutiveBurst_durations); % Mean duration


% Print calculated values
disp(['...' ]);
disp(['Burst Frequency (Hz): ' num2str(burstparameters.meanburstfrequency)]);
disp(['Interburst interval (sec): ' num2str(burstparameters.meanibi)]);
disp(['Interburst interval coefficent of variation: ' num2str(burstparameters.ibi_CV)]);
disp(['Burst peak size (spike number): ' num2str(burstparameters.meanburstsize)]);
disp(['Burst size coefficient of variation: ' num2str(burstparameters.burstsize_cv)]);
disp(['Number of Consecutive Bursts: ' num2str(burstparameters.mean_consecutiveBurst_lengths)]);
disp(['SAT Duration (sec): ' num2str(burstparameters.mean_consecutiveBurst_durations)]);



%% ============================================================
%  Visualization of detected bursts and detected parameters
%  ============================================================

figure('position', [200 300 1500 600]); hold on;
h1 = plot(time, N2(1:bin_number), 'k', 'linewidth', 1); % Plot smoothed population activity
h2 = plot([0 total_duration], [burstthreshold burstthreshold], 'r--', 'linewidth',1.5); % Plot burst detection threshold
yl = ylim;  % Get y-axis limits for shading
for i = 1:numel(runStarts) % Shade consecutive burst (SAT) periods
    x1 = burstparameters.burstlocation(runStarts(i));
    x2 = burstparameters.burstlocation(runEnds(i) + 1);
    h_temp = patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.9 1],'FaceAlpha', 0.3, 'EdgeColor', 'none');
    if i == 1
        h3 = h_temp; % Only first patch for legend
    end
end
plot(time, N2(1:bin_number), 'k', 'linewidth', 1); % Re-plot signal on top of patches
h4 = plot(burstparameters.burstlocation, burst_spikes, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); % Plot burst peaks
xlim([0 total_duration/2]); xlabel('Time (sec)'); ylabel('Spike count');
title('Detected Bursts and Consecutive Burst Periods');
set(gca, 'fontsize', 14, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5, 'tickdir', 'out');
legend([h1, h2, h3, h4],{'Smoothed Activity','Burst Threshold','Consecutive Bursts (SAT)','Burst Peaks'},'Location','northeast');% Add legend

   
   