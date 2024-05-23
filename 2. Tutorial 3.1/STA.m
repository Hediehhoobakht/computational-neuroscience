function [sta, tcorr] = STA(Iapp, spikes, dt, tminus, tplus)
% Computes the spike-triggered average.    

if (~exist('tminus'))
    tminus = 75e-3;
end
if (~exist('tplus'))
    tplus = 25e-3;
end
nminus = ceil(tminus/dt); % Number of time points before zero
nplus = ceil(tplus/dt);   % Number of time points after zero
nt = length(Iapp);     % length of original data set
sum_I = zeros(1,nminus+nplus+1);  % STA will accumulate here
tcorr = -nminus*dt:dt:nplus*dt;   % Vector of time points for STA
Iapp = Iapp - mean(Iapp); % Removes mean applied current
spikeposition = find(spikes);  % Time bins for each spike
totalspikes = length(spikeposition);    % Total number of spikes
for spike = 1:totalspikes
    ispike = spikeposition(spike);     % ispike is the bin containing a spike
    imin = max(1,ispike-nminus);       % Bin to start measuring stimulus
    imax = min(nt,ispike+nplus);       % Bin to finish measuring
    % The following lines put the stimulus, Iapp, into bins shifted
    % by the spike time (ispike)
    for i = imin:imax
        sum_I(i-ispike+nminus+1) = sum_I(i-ispike+nminus+1) ...
            + Iapp(i)/totalspikes;
    end
end

% Finally normalize by the number of contributing spikes to obtain the
% average (the mean)
sta = sum_I/totalspikes;