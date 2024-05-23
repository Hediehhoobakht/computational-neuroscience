function spikes = findPeaks(data)
  
  spikes = zeros(size(data)); 
  for i = 2:length(data)
    if data(i) < data(i-1)
      spikes(i-1) = 1;
    end
  end
end