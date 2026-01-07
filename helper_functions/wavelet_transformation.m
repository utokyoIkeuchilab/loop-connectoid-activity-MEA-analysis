function []=wavelet_transformation(Fs, time_ms, num_electrode, LP_Signal_fix, Downsample_rate, t1, t2)

%WAVELET_TRANSFORMATION  Time–frequency analysis using continuous wavelet
%transformation
%
%   INPUTS:
%       Fs              - Sampling frequency (Hz)
%       time_ms         - Time vector (milliseconds)
%       num_electrode   - Number of electrodes
%       LP_Signal_fix   - Low-pass filtered signal (< 1000 Hz)
%       Downsample_rate - Downsampling factor
%       t1, t2          - Start and end indices (after downsampling)
%
%   OUTPUT:
%       Figure generation
%
%   Description:
%       This function performs a continuous wavelet transform (CWT)
%       on low-pass filtered signals from two electrodes and plots:
%           1) The time-domain signal
%           2) The corresponding scalogram (time–frequency power)

dFs=Fs/Downsample_rate; % Effective sampling frequency after downsampling
DS_LP_Signal= downsample(LP_Signal_fix, Downsample_rate); % Downsample signals
DS_time_ms= downsample(time_ms, Downsample_rate); % Downsample time vector
DS_time_ms_lim= DS_time_ms(t1:t2,:); % Restrict time vector to analysis window

% Wavelet transform and plotting for all electrodes
for i=1:num_electrode
    DS_LP_Signal_single=DS_LP_Signal(:, i);
    
    DS_LP_Signal_lim= DS_LP_Signal_single(t1:t2,:);
    [wt,f] = cwt(DS_LP_Signal_lim, dFs,'FrequencyLimits',[0 100]);
    
% Create figure
    fig1 = figure;
    fig1.PaperUnits      = 'centimeters';
    fig1.Units           = 'centimeters';
    fig1.Color           = 'w';
    fig1.InvertHardcopy  = 'off';
    fig1.Name            = ['#', num2str(i), ' Time-Frequency analysis'];
    fig1.DockControls    = 'on';
    fig1.WindowStyle    = 'docked';
    fig1.NumberTitle     = 'off';
    set(fig1,'defaultAxesXColor','k');
    figure(fig1);
    title('Signal and Scalogram')
    
% Plot time-domain signal    
    subplot(211) 
    plot(DS_time_ms_lim/1000, DS_LP_Signal_lim);
    xlim([t1/1000 t2/1000]);
    ylim([-0.05 0.1]);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    
% Plot scalogram (time–frequency power)    
    subplot(212)
    [minf,maxf] = cwtfreqbounds(numel(DS_LP_Signal_lim),dFs);  % Determine valid frequency bounds for plotting
    AX = gca;
    freq = 2.^(round(log2(minf)):round(log2(maxf))); % Use log-spaced frequency ticks for readability
    AX.YTickLabelMode = 'auto';
    AX.YTick = freq;
    caxis([0 0.005])
    surface(DS_time_ms_lim/1000,f,abs(wt)); % Plot wavelet magnitude
    axis tight
    shading flat
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    set(gca,'yscale','log')
    colormap jet;
    
end