close all
clear
clc
Vth=-50e-3;
Vreset=-80e-3;
sigv=1e-3;
tau=3e-3;
E_L=-70e-3;
E_I=-65e-3;
E_E=0;
G_L=50e-12;
W_EE=25e-9;
W_EI=4e-9;
W_IE=800e-9;
G_in1=1e-9;
G_in2=0;
dt=0.1e-3;
t=0:dt:2.5;
tau_E=2e-3;
tau_I=5e-3;
alpha=0.2;
G_E1(1)=0;
G_E2(1)=0;
G_I1(1)=0;
G_I2(1)=0;
r_1(1)=0;
r_2(1)=0;
s_E1(1)=0;
s_I2(1)=0;
for n = 2:length(t)
     G_E1(n) = W_EE*s_E1(n-1) + G_in1;
      G_E2(n) = W_EI*s_E1(n-1) + G_in2;
      G_I1(n) = W_IE*s_I2(n-1);
      G_I2(n) = 0;
      V_ss1(n) = (G_L*E_L + G_I1(n)*E_I + G_E1(n)*E_E)/(G_L + G_I1(n) + G_E1(n));
      V_ss2(n) = (G_L*E_L + G_I2(n)*E_I + G_E2(n)*E_E)/(G_L + G_I2(n) + G_E2(n));
      %update diff eqs vie Euler's method
      r_1(n) =  r_1(n-1)+(dt/tau)*(-r_1(n-1)+fin(V_ss1(n)));
      r_2(n) =  r_2(n-1)+(dt/tau)*(-r_2(n-1)+fin(V_ss2(n)));
      s_E1(n) =  s_E1(n-1)+dt*((-s_E1(n-1)/tau_E)+alpha*r_1(n)*(1-s_E1(n-1)));
      s_I2(n) =  s_I2(n-1)+dt*((-s_I2(n-1)/tau_I)+alpha*r_2(n)*(1-s_I2(n-1)));
end
f1=figure(1);
plot(t,r_1)

hold on
plot(t,r_2)
xlabel('time');
ylabel('Firing rate');
title('part 1');
legend('r1','r2')
saveas(f1, sprintf('1.png'));
peaks_r1 = findPeaks(r_1);
peaks_r2 = findPeaks(r_2);

spike_times1 = find(peaks_r1)*dt;
spike_times2 = find(peaks_r2)*dt;
frequency_r1 = 1/mean(diff(spike_times1));
frequency_r2 = 1/mean(diff(spike_times2));
% Display the oscillation frequencies
message = sprintf('Firing rate 1 oscillation frequency: %.2f Hz', frequency_r1);
disp(message);
%Firing rate 1 oscillation frequency: 7270.65 Hz
message = sprintf('Firing rate 1 oscillation frequency: %.2f Hz', frequency_r2);
disp(message);
%Firing rate 1 oscillation frequency: 8100.11 Hz

%% Part 3
new_r1=r_1(find(t>=0.5));
new_r2=r_2(find(t>=0.5));    
[p,F]=periodogram(new_r1,[],[],1/dt);
[p2,F2]=periodogram(new_r2,[],[],1/dt);

f2=figure(2);
plot(F(1:500),p(1:500))
xlabel('frequency')
ylabel('periodgram r')
[peak loc]=max(p)
disp(F(loc))
saveas(f2, sprintf('2.png'));
% f4=figure(4);
% plot(F2(1:500),p2(1:500))
% xlabel('frequency')
% ylabel('periodgram r')
% [peak loc]=max(p)
%The f=0 component represents the mean or DC offset of the signal and has the greatest amplitude.Even if the signal oscillates around this mean, the amplitude of the oscillations is typically smaller than the mean value itself. 
%we truncate the r1 signal to remove the non-oscilatory part of signal.multiple peaks are representing a frequency at which there is significant power in the firing rate data.The number and location of the peaks are depend on the intrinsic dynamics of the coupled oscillator system

%% part 4

% Define parameters
G_in2 = 0; % Input to inhibitory cells (constant)
G_in1_range = 0:0.1e-9:10e-9; % Range of input to excitatory cells
G_E1(1)=0;
G_E2(1)=0;
G_I1(1)=0;
G_I2(1)=0;
r_1(1)=0;
r_2(1)=0;
s_E1(1)=0;
s_I2(1)=0;
i=1;
G_in1_range= 0:0.1e-9:10e-9;
% Simulation loop
for i=1:length(G_in1_range)
    G_in1=G_in1_range(i);
    % Simulate network dynamics
    for n = 2:length(t)
        G_E1(n) = W_EE*s_E1(n-1) + G_in1;
        G_E2(n) = W_EI*s_E1(n-1) + G_in2;
        G_I1(n) = W_IE*s_I2(n-1);
        G_I2(n) = 0;
        V_ss1(n) = (G_L*E_L + G_I1(n)*E_I + G_E1(n)*E_E)./(G_L + G_I1(n) + G_E1(n));
        
        
        V_ss2(n) = (G_L*E_L + G_I2(n)*E_I + G_E2(n)*E_E)./(G_L + G_I2(n) + G_E2(n));
        %update diff eqs vie Euler's method
        
        r_1(n) =  r_1(n-1)+(dt/tau)*(-r_1(n-1)+fin(V_ss1(n)));
        r_2(n) =  r_2(n-1)+(dt/tau)*(-r_2(n-1)+fin(V_ss2(n)));
        s_E1(n) =  s_E1(n-1)+dt*((-s_E1(n-1)/tau_E)+alpha*r_1(n)*(1-s_E1(n-1)));
        s_I2(n) =  s_I2(n-1)+dt*((-s_I2(n-1)/tau_I)+alpha*r_2(n)*(1-s_I2(n-1)));
    end
    new_r1=r_1(find(t>=0.5));
    new_r2=r_2(find(t>=0.5));
    new_r1=new_r1-mean(new_r1);
    new_r2=new_r2-mean(new_r2);
    [p1,F1]=periodogram(new_r1,[],[],1/dt);
    [p2,F2]=periodogram(new_r2,[],[],1/dt);
    
    [max_power, max_index] = max(p1);

% Extract the frequency corresponding to the maximum power
    frequency_at_max_power = F1(max_index);

    oscillation_frequency1(i)=frequency_at_max_power;
    [max_power2, max_index2] = max(p2);
    frequency_at_max_power2 = F2(max_index2);
    oscillation_frequency2(i)=frequency_at_max_power2;
    peaks_r1 = max(p1)-min(p1);
    peaks_r2 =  max(p2)-min(p2);
    excitatory_amplitude1(i) = peaks_r1;
    inhibitory_amplitude2(i) = peaks_r2;
    
    
%     mean_excitatory(i) = mean(new_r1> fin(Vth));     % Mean firing rate of cell 1
%     mean_inhibitory(i) = mean(new_r2> fin(Vth));     % Mean firing rate of cell 2
    mean_excitatory(i) = mean(p1);     % Mean firing rate of cell 1
    mean_inhibitory(i) = mean(p2);

end

% Plot figures
f5=figure;

subplot(2,2,1);
plot(G_in1_range, oscillation_frequency1);
hold on
plot(G_in1_range, oscillation_frequency2);
xlabel('stimulus amplitude');
ylabel('Oscillation Frequency (Hz)');
legend('E', 'I');
title('Oscillation Frequency vs. G_{in1}');

subplot(2,2,2);
plot(G_in1_range*1e9, excitatory_amplitude1);
hold on;
plot(G_in1_range*1e9, inhibitory_amplitude2);
xlabel('stimulus amplitude');
ylabel('Oscillation Amplitude');
legend('E', 'I');
title('Oscillation Amplitude vs. G_{in1}');

subplot(2,2,3);
plot(G_in1_range*1e9, mean_excitatory);
hold on;
plot(G_in1_range*1e9, mean_inhibitory);
xlabel('stimulus amplitude');
ylabel('Mean Firing Rate');
legend('E', 'I');
title('Mean Firing Rate vs. G_{in1}');
saveas(f5, sprintf('3.png'));