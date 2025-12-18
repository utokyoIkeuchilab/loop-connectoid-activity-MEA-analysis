function [LP_Signal_fix, HP_Signal_fix]=filter_signal(Fs, num_electrode, Signal)
nm=num_electrode;
time_ms = Signal(:,1);
tl=length(time_ms);

Signal_fix=zeros(tl, nm);
HPt_Signal_fix=zeros(tl, nm);

parfor i=1:nm
    
    baseline= mean(Signal(:,i));
    Signal_fix(:, i) = Signal(:, i) -baseline;
    LP_Signal_fix(:, i)=lowpass(Signal_fix(:, i),1000, Fs);
    HPt_Signal_fix(:, i)=lowpass(Signal_fix(:, i),3000, Fs);
    HP_Signal_fix(:, i)=highpass(HPt_Signal_fix(:, i),300, Fs);
end


clearvars HPt_Signal_fix Signal_fix;
end