function []=wavelet_coherence(Fs, time_ms,LP_Signal_fix, Downsample_rate, t1, t2, e1, e2, PhaseDisplayThreshold)

%WAVELET_COHERENCE  Compute and plot wavelet coherence between two electrodes
%
%   INPUTS:
%       Fs                    - Original sampling frequency (Hz)
%       time_ms               - Time vector (milliseconds)
%       LP_Signal_fix         - Low-pass filtered signal (<1000 Hz)
%       Downsample_rate       - Downsampling factor
%       t1, t2                - Start and end indices (after downsampling)
%       e1, e2                - Electrode indices to compare
%       PhaseDisplayThreshold - Coherence threshold for phase arrows
%
%   OUTPUT:
%       Figure generation
%
%   Description:
%       This function downsamples the signals from two electrodes,
%       and computes their wavelet coherence over a selected time window.

dFs=Fs/Downsample_rate; % Effective sampling frequency after downsampling
DS_LP_Signal= downsample(LP_Signal_fix, Downsample_rate); % Downsample signals
DS_time_ms= downsample(time_ms, Downsample_rate); % Downsample time vector
DS_time_ms_lim= DS_time_ms(t1:t2,:); % Limit time vector to analysis window

% Extract signals from selected electrodes
DS_LP_Signal_single_e1=DS_LP_Signal(:, e1); 
DS_LP_Signal_single_e2=DS_LP_Signal(:, e2);
DS_LP_Signal_single_e1_lim= DS_LP_Signal_single_e1(t1:t2,:);
DS_LP_Signal_single_e2_lim= DS_LP_Signal_single_e2(t1:t2,:);

% Create figure
    fig1 = figure;
    fig1.PaperUnits      = 'centimeters';
    fig1.Units           = 'centimeters';
    fig1.Color           = 'w';
    fig1.InvertHardcopy  = 'off';
    fig1.Name            = [num2str(e1), 'vs', num2str(e2),' Wavelet_coherence'];
    fig1.DockControls    = 'on';
    fig1.WindowStyle    = 'docked';
    fig1.NumberTitle     = 'off';
    set(fig1,'defaultAxesXColor','k');
    figure(fig1);

% Compute wavelet coherence
    wcoherence(DS_LP_Signal_single_e1_lim,DS_LP_Signal_single_e2_lim,dFs, 'PhaseDisplayThreshold',PhaseDisplayThreshold);
    ax = gca;
    ytick=round(pow2(ax.YTick),3);
    ax.YTickLabel=ytick;
    ax.YLabel.String='Frequency (Hz)';
    ax.Title.String = 'Wavelet Coherence';
    colormap(jet)
    caxis([0 1.5])
end