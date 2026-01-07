function[Pos_all_spikes, Neg_all_spikes, pos_avg_amp, neg_avg_amp, num_posspks, num_negspks, All_interspike_interval_sec, interspike_interval_sec_avg, All_spikes]=spike_detection(Fs, time_ms, num_electrode, HP_Signal_fix, v, mag)

%SPIKE_DETECTION  Detect positive and negative spikes from MEA signals
%
%   INPUTS:
%       Fs              - Sampling frequency (Hz)
%       time_ms         - Time vector in milliseconds
%       num_electrode   - Number of electrodes
%       HP_Signal_fix   - High-pass filtered signal (300–3000 Hz)
%       v               - Visualization flag (1 = ON, 0 = OFF)
%       mag             - Spike detection threshold multiplier (× STD)
%
%   OUTPUTS:
%       Pos_all_spikes  - Cell array:
%                           {i,1} = positive spike times (sec)
%                           {i,2} = positive spike amplitudes
%       Neg_all_spikes  - Cell array:
%                           {i,1} = negative spike times (sec)
%                           {i,2} = negative spike amplitudes
%       pos_avg_amp     - Mean positive spike amplitude per electrode
%       neg_avg_amp     - Mean negative spike amplitude per electrode
%       num_posspks     - Number of positive spikes per electrode
%       num_negspks     - Number of negative spikes per electrode
%       All_interspike_interval_sec
%                       - Cell array of ISIs (sec) per electrode
%       interspike_interval_sec_avg
%                       - Mean ISI per electrode (sec)
%       All_spikes      - Cell array:
%                           {i,1} = all spike times (sec)
%                           {i,2} = electrode index for raster plots
%   
%   Detection method:
%       - Threshold = mag × standard deviation of signal
%       - Positive and negative peaks detected separately
%       - Interspike intervals computed from combined spikes

if v==1
 fig1 = figure;
 fig1.PaperUnits      = 'centimeters';
 fig1.Units           = 'centimeters';
 fig1.Color           = 'w';
 fig1.InvertHardcopy  = 'off';
 fig1.Name            = 'Spike detection overview';
 fig1.NumberTitle     = 'off';
 set(fig1,'defaultAxesXColor','k');
 figure(fig1);
end

Pos_all_spikes              = {};
Neg_all_spikes              = {};
All_spikes                  = {};
pos_avg_amp                 = zeros(num_electrode, 1);
neg_avg_amp                 = zeros(num_electrode, 1);
interspike_interval_sec_avg = zeros(num_electrode, 1);
num_posspks                 = zeros(num_electrode, 1);
num_negspks                 = zeros(num_electrode, 1);

for i=1:num_electrode
    STDEV=std(HP_Signal_fix(:,i));
    peak_th=mag*STDEV; % set threshold
    [posspks, poslocs] = findpeaks(HP_Signal_fix(:,i), Fs,'MinPeakHeight',peak_th ); % Detect positive spikes
    [negspks, neglocs] = findpeaks(-HP_Signal_fix(:,i), Fs,'MinPeakHeight',peak_th ); % Detect negative spikes
    
    all_locs=vertcat(poslocs,neglocs);  % spike times (sec)
    all_spks=vertcat(posspks,-negspks);  % signed amplitudes
    
    % Sort spikes by time
    temp_all_locs_spks=horzcat(all_locs, all_spks);
    all_locs_spks=sortrows(temp_all_locs_spks, 1);
    interspike_interval=diff(all_locs_spks(:,1)); % Interspike intervals (sec)
    
    Pos_all_spikes{i,1}=poslocs;
    Neg_all_spikes{i,1}=neglocs;
    Pos_all_spikes{i,2}=posspks;
    Neg_all_spikes{i,2}=negspks;
    
    temp_locs=sort(vertcat(poslocs, neglocs));
    num_locs=length(temp_locs);
    dummy_mat= ones(num_locs,1);
    All_spikes{i,1}=temp_locs;
    All_spikes{i,2}=dummy_mat*i;
    
    All_interspike_interval_sec{i,1}=interspike_interval;
    interspike_interval_sec_avg(i,1)=mean(interspike_interval);
    
    pos_avg_amp(i, 1)=mean(posspks);
    neg_avg_amp(i, 1)=mean(-negspks);
    num_posspks(i, 1)=length(posspks);
    num_negspks(i, 1)=length(negspks);
    
    % Optional visualization
    if v==1
        subplot(num_electrode, 1, i)
        hold on
        plot(time_ms/1000, HP_Signal_fix(:, i), 'Color', 'black');
        neglocs=Neg_all_spikes{i,1};
        negspks=Neg_all_spikes{i,2};
        poslocs=Pos_all_spikes{i,1};
        posspks=Pos_all_spikes{i,2};
        plot(neglocs, -negspks, 'o');
        plot(poslocs, posspks, 'o');
        yline(peak_th,'r');
        yline(-peak_th,'r');
        xlim([0 30]);
        hold off
        ylabel('Amplitude (mV)')
        xlabel('Time (sec)')
     end


end
clearvars temp_all_locs_spks all_locs  all_spks all_locs_spks
