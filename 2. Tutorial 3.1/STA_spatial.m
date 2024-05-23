function [sta, tcorr] = STA_spatial(stim_array, spikes, dt,new_nstep_length, tminus, tplus)
  % Computes the spatiotemporal spike-triggered average.

  % Handle default values
  if (~exist('tminus'))
    tminus = 75e-3;
  end
  if (~exist('tplus'))
    tplus = 25e-3;
  end

  % Get dimensions and calculate number of time points before/after zero
  [Nspace, Nt] = size(stim_array);
  nminus = ceil(tminus/dt);
  nplus = ceil(tplus/dt);

  % Initialize STA and time vector
  sta = zeros(Nspace, nminus+nplus+1);
  tcorr = -nminus*dt:dt:nplus*dt;

  % Remove mean from each spatial bin
  for step = 1:Nt/new_nstep_length;                    % Loop through steps of constant current
      step
    istart = (step-1)*new_nstep_length+1; % first time point with 1ms steps
    istop = step*new_nstep_length;        % last time point with 1ms steps
    istart
    istop
    % generate the input vector as before, but with bins of 1ms, not of dt
    stim_array(:,istart:istop) =  repmat((stim_array(:,istart)- mean(stim_array, 2)),1,new_nstep_length);
end
%   for i=1:Nt
%    stim_array(:,i) = stim_array(:,i) - mean(stim_array, 2);
%   end
  % Find spike positions and total number of spikes
  spikeposition = find(spikes);
  totalspikes = length(spikeposition);

  % Loop through spikes
  for spike = 1:totalspikes
    ispike = spikeposition(spike); % Time bin containing the spike
    imin = max(1, ispike-nminus); % Start of window
    imax = min(Nt, ispike+nplus); % End of window

    % Accumulate stimulus values for each spatial bin
    for i = imin:imax
      for j = 1:Nspace
        sta(j, i-ispike+nminus+1) = sta(j, i-ispike+nminus+1) + stim_array(j, i)/totalspikes;
      end
    end
  end

  % Normalize by the number of contributing spikes
  sta = sta / totalspikes;

  % Return STA and time vector
  return;
end
