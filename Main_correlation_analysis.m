%% ============================================================
%  Main correlation analysis script
%  - Aggregate spike trains across electrodes
%  - Compute population firing-rate time series
%  - Estimate inter-organoid and intra-organoid correlations
%  ============================================================

clear all
close all
clc



%% ============================================================
%  Set parameters and import data
%  ============================================================

total_duration = 1200; % Total recording duration (seconds)
bin_win = 100; % Binning window (ms) for spike counts
smoothingfactor = 3;  % Gaussian smoothing strength
bin_number = total_duration*1000/bin_win; % Number of bins for the full recording
time = [total_duration/bin_number:total_duration/bin_number:total_duration];

% Load MATLAB .mat file containing the spike timing data for all 64 electrodes (can be downloaded on Github)
All_spikes = struct2array(load('C:\Users\Tomoya\Projects\papers\loop-connectoid-activity-MEA-analysis\example_data\example_spike_data.mat'));


% Define electrode location
TL = [1,2,3,4,9,10,11,12,17,18,19,20,25,26,27,28]; % Electrodes of top left organoids
TR = [5,6,7,8,13,14,15,16,21,22,23,24,29,30,31,32]; % Electrodes of top right organoids
BL = [33,34,35,36,41,42,43,44,49,50,51,52,57,58,59,60]; % Electrodes of bottom left organoids
BR = [37,38,39,40,45,46,47,48,53,54,55,56,61,62,63,64]; % Electrodes of bottom right organoids



%% ============================================================
%  Population spike aggregation and binning
%  ============================================================

% Concatenate all spike times across all 16 electrodes per organoid
SpikesTL = [];
SpikesTR = [];
SpikesBL = [];
SpikesBR = [];
for i = 1:16 % write all spike timings into one vector
    SpikesTL = [SpikesTL; All_spikes{TL(i), 1}];
    SpikesTR = [SpikesTR; All_spikes{TR(i), 1}];
    SpikesBL = [SpikesBL; All_spikes{BL(i), 1}];
    SpikesBR = [SpikesBR; All_spikes{BR(i), 1}];
end

% Bin spike counts over time
[NTL,~] = histcounts(SpikesTL, bin_number);
[NTR,~] = histcounts(SpikesTR, bin_number);
[NBL,~] = histcounts(SpikesBL, bin_number);
[NBR,~] = histcounts(SpikesBR, bin_number);

% Smooth binned spike counts using a Gaussian kernel
NTL = conv(NTL,gausswin(smoothingfactor*(1/(bin_win/1000))));
NTR = conv(NTR,gausswin(smoothingfactor*(1/(bin_win/1000))));
NBL = conv(NBL,gausswin(smoothingfactor*(1/(bin_win/1000))));
NBR = conv(NBR,gausswin(smoothingfactor*(1/(bin_win/1000))));



%% ============================================================
%  Inter-organoid correlation for 4 connected organoids
%  ============================================================

% Correlation between population activity of the four organoids
correlation = corrcoef([NTL' NTR' NBL' NBR']);

% Selected correlations (direct connections between the 4 organoids)
intercorrelation = [correlation(1,2); correlation(1,3); correlation(2,4); correlation(3,4)];
disp('Inter organoid correlations')
disp(['Inter organoid correlation (direct connections): ' num2str(correlation(1,2))])
disp(['Inter organoid correlation (direct connections): ' num2str(correlation(1,3))])
disp(['Inter organoid correlation (direct connections): ' num2str(correlation(2,4))])
disp(['Inter organoid correlation (direct connections): ' num2str(correlation(3,4))])

% Diagonal (indirect) organoid correlations
indirectcorrelation = [correlation(1,4); correlation(2,3);];
disp(['Inter organoid correlation (indirect connections): ' num2str(correlation(1,4))])
disp(['Inter organoid correlation (indirect connections): ' num2str(correlation(2,3))])

% Plot population activity traces of the 4 organoids
figure('position', [300 300 1200 600]); hold on;
plot(time,NTL(1:bin_number), 'color', 'k', 'linewidth',1.5); plot(time,NTR(1:bin_number), 'color', 'b', 'linewidth',1.5);
plot(time,NBL(1:bin_number), 'color', 'r', 'linewidth',1.5); plot(time,NBR(1:bin_number), 'color', 'g', 'linewidth',1.5);
xlim([0 total_duration]); ylabel('Spike count'); xlabel('Time (sec)'); title('Inter organoid activity');
set(gca, 'fontsize', 14, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5);
legend({'Top Left Organoid' 'Top Right Organoid' 'Bottom Left Organoid' 'Bottom Right Organoid'}, 'box', 'off')

% Plot correlation matrix
figure
imagesc(correlation); c = colorbar;
c.Ticks = [0 0.2 0.4 0.6 0.8 1];
c.TickLabels = {'0' '' '' '' '' '1.0' }; caxis([0 1]);
c.TickDirection = 'out'; c.LineWidth = 1.5;
ylabel(c,'Correlation Coefficient','fontsize', 20, 'fontname','san serif', 'fontweight', 'bold','Rotation',270);
set (gca, 'fontsize', 22, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5);
set(gca,'xtick',[], 'ytick', []);
ylabel('Organoid #'); xlabel('Organoid #'); title('Inter Organoid Correlation','fontsize', 20)



%% ============================================================
%  Intra-organoid correlation (within-organoid)
%  ============================================================

intraelectrodenumbers = 3; % Number of highest-rate electrodes to include

% Compute firing rates per electrode for each organoid
elespikerateTL = [];
elespikerateTR = [];
elespikerateBL = [];
elespikerateBR = [];
for i = 1:16
    elespikerateTL = [elespikerateTL length(All_spikes{TL(i), 1})/total_duration];
    elespikerateTR = [elespikerateTR length(All_spikes{TR(i), 1})/total_duration];
    elespikerateBL = [elespikerateBL length(All_spikes{BL(i), 1})/total_duration];
    elespikerateBR = [elespikerateBR length(All_spikes{BR(i), 1})/total_duration];
end

%  Extract top-N (3) electrodes and compute spiking activity
for i = 1:intraelectrodenumbers
    [~,idx] = maxk(elespikerateTL,intraelectrodenumbers);
    [N,~] = histcounts([All_spikes{TL(idx(i)), 1};0;total_duration], bin_number);
    N = conv(N,gausswin(smoothingfactor*(1/(bin_win/1000))));
    IntraSpikesTL(i,:) = N;

    [~,idx] = maxk(elespikerateTR,intraelectrodenumbers);
    [N,~] = histcounts([All_spikes{TR(idx(i)), 1};0;total_duration], bin_number);
    N = conv(N,gausswin(smoothingfactor*(1/(bin_win/1000))));
    IntraSpikesTR(i,:) = N;
    
    [~,idx] = maxk(elespikerateBL,intraelectrodenumbers);
    [N,~] = histcounts([All_spikes{BL(idx(i)), 1};0;total_duration], bin_number);
    N = conv(N,gausswin(smoothingfactor*(1/(bin_win/1000))));
    IntraSpikesBL(i,:) = N;

    [~,idx] = maxk(elespikerateBR,intraelectrodenumbers);
    [N,~] = histcounts([All_spikes{BR(idx(i)), 1};0;total_duration], bin_number);
    N = conv(N,gausswin(smoothingfactor*(1/(bin_win/1000))));
    IntraSpikesBR(i,:) = N;

end



%% ============================================================
%  Intra-organoid correlation plots and averages
%  ============================================================
disp('...')
disp('Intra organoid correlations')


% --- Top left organoid ---
correlation = corrcoef(IntraSpikesTL');
mean_intra_corr1 = mean([correlation(2,1); correlation(3,1); correlation(3,2)]);
disp(['Intra correlation top left organoid: ' num2str(mean_intra_corr1)])

% Plot spiking activity of electrodes within organoid 
figure('position', [300 300 1200 600]); hold on;
plot(time,IntraSpikesTL(:,1:bin_number), 'linewidth',1.5)
xlim([0 total_duration]); ylabel('Spike count'); xlabel('Time (sec)'); title('Intra organoid activity top left');
set(gca, 'fontsize', 14, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5);

figure; imagesc(correlation); c = colorbar;
c.Ticks = [0 0.2 0.4 0.6 0.8 1];
c.TickLabels = {'0' '' '' '' '' '1.0' }; caxis([0 1]);
c.TickDirection = 'out'; c.LineWidth = 1.5;
ylabel(c,'Correlation Coefficient','fontsize', 20, 'fontname','san serif', 'fontweight', 'bold','Rotation',270);
set (gca, 'fontsize', 22, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5);
set(gca,'xtick',[], 'ytick', []);
ylabel('Electrode #'); xlabel('Electrode #'); title('Intra Organoid Correlation Top Left','fontsize', 16)



% --- Top right organoid ---
correlation = corrcoef(IntraSpikesTR');
mean_intra_corr2 = mean([correlation(2,1); correlation(3,1); correlation(3,2)]);
disp(['Intra correlation top right organoid: ' num2str(mean_intra_corr2)])

% Plot spiking activity of electrodes within organoid 
figure('position', [300 300 1200 600]); hold on;
plot(time,IntraSpikesTR(:,1:bin_number), 'linewidth',1.5)
xlim([0 total_duration]); ylabel('Spike count'); xlabel('Time (sec)'); title('Intra organoid activity top right');
set(gca, 'fontsize', 14, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5);

% Plot correlation matrix
figure; imagesc(correlation); c = colorbar;
c.Ticks = [0 0.2 0.4 0.6 0.8 1];
c.TickLabels = {'0' '' '' '' '' '1.0' }; caxis([0 1]);
c.TickDirection = 'out'; c.LineWidth = 1.5;
ylabel(c,'Correlation Coefficient','fontsize', 20, 'fontname','san serif', 'fontweight', 'bold','Rotation',270);
set (gca, 'fontsize', 22, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5);
set(gca,'xtick',[], 'ytick', []);
ylabel('Electrode #'); xlabel('Electrode #'); title('Intra Organoid Correlation Top Right','fontsize', 16)



% --- Bottom left organoid ---
correlation = corrcoef(IntraSpikesBL');
mean_intra_corr3 = mean([correlation(2,1); correlation(3,1); correlation(3,2)]);
disp(['Intra correlation bottom left organoid: ' num2str(mean_intra_corr3)])

% Plot spiking activity of electrodes within organoid 
figure('position', [300 300 1200 600]); hold on;
plot(time,IntraSpikesBL(:,1:bin_number), 'linewidth',1.5)
xlim([0 total_duration]); ylabel('Spike count'); xlabel('Time (sec)'); title('Intra organoid activity bottom left');
set(gca, 'fontsize', 14, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5);

% Plot correlation matrix
figure; imagesc(correlation); c = colorbar;
c.Ticks = [0 0.2 0.4 0.6 0.8 1];
c.TickLabels = {'0' '' '' '' '' '1.0' }; caxis([0 1]);
c.TickDirection = 'out'; c.LineWidth = 1.5;
ylabel(c,'Correlation Coefficient','fontsize', 20, 'fontname','san serif', 'fontweight', 'bold','Rotation',270);
set (gca, 'fontsize', 22, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5);
set(gca,'xtick',[], 'ytick', []);
ylabel('Electrode #'); xlabel('Electrode #'); title('Intra Organoid Correlation Bottom Left','fontsize', 16)



% --- Bottom right organoid ---
correlation = corrcoef(IntraSpikesBR');
mean_intra_corr4 = mean([correlation(2,1); correlation(3,1); correlation(3,2)]);
disp(['Intra correlation bottom right organoid: ' num2str(mean_intra_corr4)])

% Plot spiking activity of electrodes within organoid 
figure('position', [300 300 1200 600]); hold on;
plot(time,IntraSpikesBR(:,1:bin_number), 'linewidth',1.5)
xlim([0 total_duration]); ylabel('Spike count'); xlabel('Time (sec)'); title('Intra organoid activity bottom right');
set(gca, 'fontsize', 14, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5);

% Plot correlation matrix
figure; imagesc(correlation); c = colorbar;
c.Ticks = [0 0.2 0.4 0.6 0.8 1];
c.TickLabels = {'0' '' '' '' '' '1.0' }; caxis([0 1]);
c.TickDirection = 'out'; c.LineWidth = 1.5;
ylabel(c,'Correlation Coefficient','fontsize', 20, 'fontname','san serif', 'fontweight', 'bold','Rotation',270);
set (gca, 'fontsize', 22, 'fontname','san serif', 'fontweight', 'bold', 'linewidth', 1.5);
set(gca,'xtick',[], 'ytick', []);
ylabel('Electrode #'); xlabel('Electrode #'); title('Intra Organoid Correlation Bottom Right','fontsize', 16)




