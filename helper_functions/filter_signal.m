function [LP_Signal_fix, HP_Signal_fix]=filter_signal(Fs, num_electrode, Signal)

%Filtering  MEA signals
%   INPUTS:
%       Fs            - Sampling frequency (Hz)
%       num_electrode - Number of electrodes (channels)
%       Signal        - Raw signal matrix [time x electrodes]
%
%   OUTPUTS:
%       LP_Signal_fix - Low-pass filtered signal (< 1000 Hz)
%       HP_Signal_fix - Band-pass filtered signal (300–3000 Hz)
%
%   Processing steps:
%       1) Offset correction
%       2) Low-pass filtering for LFP-like signals
%       3) Band-pass filtering for spike detection

nm=num_electrode; % Number of electrodes
time_ms = Signal(:,1);
tl=length(time_ms); % Number of time samples

Signal_fix=zeros(tl, nm);
HPt_Signal_fix=zeros(tl, nm); % Temporary array for intermediate filtering

parfor i=1:nm
    
    baseline= mean(Signal(:,i));
    Signal_fix(:, i) = Signal(:, i) -baseline; % Offset correction
    LP_Signal_fix(:, i)=lowpass(Signal_fix(:, i),1000, Fs); % Low-pass filtering (< 1000 Hz)
    HPt_Signal_fix(:, i)=lowpass(Signal_fix(:, i),3000, Fs);
    HP_Signal_fix(:, i)=highpass(HPt_Signal_fix(:, i),300, Fs); % Band-pass filtering (300–3000 Hz)
end


clearvars HPt_Signal_fix Signal_fix;
end