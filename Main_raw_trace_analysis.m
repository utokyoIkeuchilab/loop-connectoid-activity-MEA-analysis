%% Main script
clear all
close all
clc

%% Set parameters
Fs=20000; % Sampling frequency

%% Load and import data
signal = struct2array(load('C:\Users\Tomoya\Projects\papers\loop-connectoid-activity-MEA-analysis\example_data\example_signal_trace.mat'));


%% Preprocess data 
num_electrode=4;

% filter Signal, LP_Signal_fix = lowpass <1000Hz, HP_Signal_fix = bandpass 300-3000Hz
[LP_Signal_fix, HP_Signal_fix]=filter_signal(Fs, num_electrode, signal);

time = [1/Fs:1/Fs:size(signal,1)/Fs]';
time_ms = time*1000;
tl=length(time_ms);


%% Spike detection
%if you need visuaiztion of spike detection result set visual_on=1
%otherwise visual_on=0

visual_on=1;
magnification=4; % magnification *STDEV for spike detection

[All_spikes_pos, All_spikes_neg, Mean_posspks_amp, Mean_negspks_amp, Num_posspks,Num_negspks, All_interspike_interval_sec, Mean_interspike_interval_sec, All_spikes]=spike_detection(Fs, time_ms, num_electrode, HP_Signal_fix, visual_on, magnification);

%% Plot all traces together

color1 = [50, 50, 180]/255;
color2 = [80, 80, 180]/255;
color3 = [110, 110, 180]/255;
color4 = [140, 140, 180]/255;

figure('position', [200 100 1500 800]) 
hold on;
h1 = plot(time, LP_Signal_fix(:,1), 'color', color1);
h2 = plot(time,LP_Signal_fix(:,2)-0.02,'color', color2);
h3 = plot(time,LP_Signal_fix(:,3)-0.03,'color', color3);
h4 = plot(time,LP_Signal_fix(:,4)-0.04,'color', color4);

xlabel('Time (s)'); ylabel('Amplitude (mV)');
set (gca, 'fontsize', 14, 'fontname','san serif', 'fontweight', 'bold', 'box', 'off', 'linewidth',1.5);
xlim([50 150])

legend([h1 h2 h3 h4], 'Organoid 1', 'Organoid 2', 'Organoid 3', 'Organoid 4', 'box', 'off')
hold off



%% Wavelet transformation of all electrodes
Downsample_rate=20;
t1 = 50; %start time of calculation range (in seconds)
t2 = 150; %end time of calculation range (in seconds)

wavelet_transformation(Fs, time_ms, num_electrode, LP_Signal_fix, Downsample_rate, t1*Fs/Downsample_rate, t2*Fs/Downsample_rate);


%% Wavelet coherence organoid 1 vs organoid 2

Downsample_rate=20;
dFs = Fs/Downsample_rate;
t1 = 1;  %start time of calculation range (in seconds)
t2 = 301; %end time of calculation range (in seconds)
e1=1; %first electrode number
e2=2; %second electrode number
PhaseDisplayThreshold=1;
DS_time_ms= downsample(time_ms(t1*Fs:t2*Fs), Downsample_rate);


wavelet_coherence(Fs, time_ms,LP_Signal_fix, Downsample_rate, t1*Fs/Downsample_rate, t2*Fs/Downsample_rate, e1, e2, PhaseDisplayThreshold)

%% Mean wavelet coherence between all 4 organoids
DS_LP_Signal= downsample(LP_Signal_fix, Downsample_rate);

% wavelet coherence between electrode 1 & 2
w1 = wcoherence(DS_LP_Signal(t1*dFs:t2*dFs,1),DS_LP_Signal(t1*dFs:t2*dFs,2),dFs, 'PhaseDisplayThreshold',PhaseDisplayThreshold);
% wavelet coherence between electrode 1 & 3
w2 = wcoherence(DS_LP_Signal(t1*dFs:t2*dFs,1),DS_LP_Signal(t1*dFs:t2*dFs,3),dFs, 'PhaseDisplayThreshold',PhaseDisplayThreshold);
% wavelet coherence between electrode 2 & 4
w3 = wcoherence(DS_LP_Signal(t1*dFs:t2*dFs,2),DS_LP_Signal(t1*dFs:t2*dFs,4),dFs, 'PhaseDisplayThreshold',PhaseDisplayThreshold);
% wavelet coherence between electrode 3 & 4
w4 = wcoherence(DS_LP_Signal(t1*dFs:t2*dFs,3),DS_LP_Signal(t1*dFs:t2*dFs,4),dFs, 'PhaseDisplayThreshold',PhaseDisplayThreshold);


figure('position', [500 500 1000 600])
imagesc((w1+w2+w3+w4)/4);
colormap(jet)
caxis([0 1]);
yticks([8:12:188])
yticklabels({'256' '128' '64' '32' '16' '8' '4' '2' '1' '0.5' '0.25' '0.125' '0.063' '0.031' '0.016' '0.008'})
xticks([1:(length(DS_LP_Signal(t1*dFs:t2*dFs,4))/2)-1:length(DS_LP_Signal(t1*dFs:t2*dFs,4))])
xticklabels({num2str(t1) num2str(round((length(DS_LP_Signal(t1*dFs:t2*dFs,4))/2)-1)/dFs) num2str(t2)})
xlabel('Time (sec)');
ylabel('Frequency (Hz)')
set (gca, 'fontsize', 16, 'fontname','san serif', 'fontweight', 'bold', 'box', 'off', 'linewidth',1.5);
set(gca,'TickLength',[0.01, 0.01])
hold off


%% Plot Different frequencies

selel = 1;

[wt,f] = cwt(DS_LP_Signal(t1*dFs:t2*dFs,selel),dFs);
xrec1 = icwt(wt,f,[0.2 0.5],'SignalMean',mean(DS_LP_Signal(t1*dFs:t2*dFs,selel)));
xrec2 = icwt(wt,f,[0.5 4],'SignalMean',mean(DS_LP_Signal(t1*dFs:t2*dFs,selel))); %delta
xrec3 = icwt(wt,f,[4 8],'SignalMean',mean(DS_LP_Signal(t1*dFs:t2*dFs,selel))); %theta
xrec4 = icwt(wt,f,[8 12],'SignalMean',mean(DS_LP_Signal(t1*dFs:t2*dFs,selel))); %alpha
xrec5 = icwt(wt,f,[12 30],'SignalMean',mean(DS_LP_Signal(t1*dFs:t2*dFs,selel))); %beta
xrec6 = icwt(wt,f,[30 300],'SignalMean',mean(DS_LP_Signal(t1*dFs:t2*dFs,selel)));
  
%Offset Correction 
xrec1 = xrec1-mean(xrec1);
xrec2 = xrec2-mean(xrec2);
xrec3 = xrec3-mean(xrec3);
xrec4 = xrec4-mean(xrec4);
xrec5 = xrec5-mean(xrec5);
xrec6 = xrec6-mean(xrec6);

figure('position', [200 100 1500 800]) 
h1 = plot(time_ms(t1*Fs:t2*Fs)/1000, LP_Signal_fix(t1*Fs:t2*Fs, selel), 'Color', [0.75 0.75 0.75]);
hold on
h6 = plot(DS_time_ms/1000, xrec6, 'Color', [1 0 0], 'linewidth',1.5)
h5 = plot(DS_time_ms/1000, xrec5, 'Color', [0 1 1], 'linewidth',1.5)
h4 = plot(DS_time_ms/1000, xrec4, 'Color', [1 1 0], 'linewidth',1.5)
h3 = plot(DS_time_ms/1000, xrec3, 'Color', [0 0 1], 'linewidth',1.5)
h2 = plot(DS_time_ms/1000, xrec2, 'Color', [0 0 0], 'linewidth',1.5)

set(gca, 'tickdir', 'out')
set (gca, 'fontsize', 14, 'fontname','san serif', 'fontweight', 'bold', 'box', 'off', 'linewidth',1.5);
xlim([50 150])
ylim([-0.05 0.07])
xlabel('Time (sec)');
ylabel('Amplitude (mV)')
legend([h1 h2 h3 h4 h5 h6], {'Raw' 'Delta' 'Theta' 'Alpha' 'Beta' 'Gamma'}, 'box', 'off')

