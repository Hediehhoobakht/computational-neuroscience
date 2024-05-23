close all
clear
clc

G_Leak = 30E-9; % Leak conductance (S)
R_m = 1/G_Leak;% resistance (ohm)
G_Na = 12E-6; %sodium conductance 
G_K = 3.6E-6; %potassium conductance
E_Na = 45E-3;
E_K = -82E-3;
leak_potential = -60E-3;
C_m = 100E-12;
tau_m = R_m*C_m;
dt = 0.02E-4;
tmax = 0.35;
time_vector =0:dt:tmax;

v = zeros(size(time_vector));
v(1) = leak_potential;
%part a
applied_current = zeros(size(time_vector));
%initialize gating variables
m = zeros(size(time_vector));
h = zeros(size(time_vector));
n = zeros(size(time_vector));
m(1) = 0;
h(1) = 0;
n(1) = 0;
for j = 1:length(time_vector)-1
    v(j+1) = v(j) + dt*( G_Leak*(leak_potential-v(j)) + G_Na*(m(j)^3)*h(j)*(E_Na-v(j))+G_K*(n(j)^4)*(E_K-v(j))+ applied_current(j))/C_m;
    if ( v(j) == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m=(1e5*(-v(j)-0.045))/(exp(100*(-v(j)-0.045))-1);
    end
    
    beta_m=4e3*exp((-v(j)-0.070)/0.018);
    alpha_h=70*exp(50*(-v(j)-0.070));
    beta_h=1e3/(1+exp(100*(-v(j)-0.040)));
    if ( v(j) == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else;                   % potassium activation rate
        alpha_n=(1e4*(-v(j)-0.060))/(exp(100*(-v(j)-0.060))-1);
    end
    
    beta_n=125*exp((-v(j)-0.070)/0.08);
    m(j+1)=m(j)+dt*(alpha_m*(1-m(j))-beta_m*m(j) );
    h(j+1)=h(j)+dt*(alpha_h*(1-h(j))-beta_h*h(j));
    n(j+1)=n(j)+dt*(alpha_n*(1-n(j))-beta_n*n(j));
end
f1=figure;
plot(time_vector, v);
xlabel('Time (seconds)');
ylabel('Voltage membrane (V)');
title('Voltage membrane vs. Time- 1.a');
grid on
saveas(f1, sprintf('1_a.png'));
%% part b
% Create applied current vector with step
v = zeros(size(time_vector));
v(1) = leak_potential;
m = zeros(size(time_vector));
h = zeros(size(time_vector));
n = zeros(size(time_vector));
m(1) = 0;
h(1) = 0;
n(1) = 0;
applied_current(time_vector >= 100e-3 & time_vector < 200e-3) = 0.22e-9;
for j = 1:length(time_vector)-1
    v(j+1) = v(j) + dt*( G_Leak*(leak_potential-v(j)) + G_Na*(m(j)^3)*h(j)*(E_Na-v(j))+G_K*(n(j)^4)*(E_K-v(j))+ applied_current(j))/C_m;
    if ( v(j) == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m=(1e5*(-v(j)-0.045))/(exp(100*(-v(j)-0.045))-1);
    end
    beta_m=4e3*exp((-v(j)-0.070)/0.018);
    alpha_h=70*exp(50*(-v(j)-0.070));
    beta_h=1e3/(1+exp(100*(-v(j)-0.040)));
    if ( v(j) == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else;                   % potassium activation rate
        alpha_n=(1e4*(-v(j)-0.060))/(exp(100*(-v(j)-0.060))-1);
    end
    beta_n=125*exp((-v(j)-0.070)/0.08);
    m(j+1)=m(j)+dt*(alpha_m*(1-m(j))-beta_m*m(j) );
    h(j+1)=h(j)+dt*(alpha_h*(1-h(j))-beta_h*h(j));
    n(j+1)=n(j)+dt*(alpha_n*(1-n(j))-beta_n*n(j));
end
% Plot the applied current
f2=figure;
subplot(2,1,1);
plot(time_vector*1e3, applied_current*1e9)
xlabel('Time (ms)')
ylabel('Applied Current (nA)')
title('Applied Current with Step - 1b')
grid on
subplot(2,1,2);
plot(time_vector*1e3, v*1e3)
xlabel('Time (ms)')
ylabel('membrane potential (mV)')
title('voltage - 1b')
grid on
saveas(f2, sprintf('1_b.png'));
%% Part c
applied_current = zeros(size(time_vector));
duration=5e-3/dt;
delay=5e-3/dt;
for j=1:10
    startpoint=((100e-3/dt)+delay*(j-1));
    endpoint=(startpoint+duration);
    applied_current(startpoint: endpoint) = 0.22e-9;
    
   
end
v = zeros(size(time_vector));
v(1) = leak_potential;
m = zeros(size(time_vector));
h = zeros(size(time_vector));
n = zeros(size(time_vector));
m(1) = 0;
h(1) = 0;
n(1) = 0;
for j = 1:length(time_vector)-1
    v(j+1) = v(j) + dt*( G_Leak*(leak_potential-v(j)) + G_Na*(m(j)^3)*h(j)*(E_Na-v(j))+G_K*(n(j)^4)*(E_K-v(j))+ applied_current(j))/C_m;
    if ( v(j) == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m=(1e5*(-v(j)-0.045))/(exp(100*(-v(j)-0.045))-1);
    end
    beta_m=4e3*exp((-v(j)-0.070)/0.018);
    alpha_h=70*exp(50*(-v(j)-0.070));
    beta_h=1e3/(1+exp(100*(-v(j)-0.040)));
    if ( v(j) == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else;                   % potassium activation rate
        alpha_n=(1e4*(-v(j)-0.060))/(exp(100*(-v(j)-0.060))-1);
    end
    beta_n=125*exp((-v(j)-0.070)/0.08);
    m(j+1)=m(j)+dt*(alpha_m*(1-m(j))-beta_m*m(j) );
    h(j+1)=h(j)+dt*(alpha_h*(1-h(j))-beta_h*h(j));
    n(j+1)=n(j)+dt*(alpha_n*(1-n(j))-beta_n*n(j));
end
f3=figure;
subplot(2,1,1);
plot(time_vector*1e3, applied_current*1e9)
xlabel('Time (ms)')
ylabel('Applied Current (nA)')
title('Applied Current with Step - 1c - 5ms delay')
grid on
subplot(2,1,2);
plot(time_vector*1e3, v*1e3)
xlabel('Time (ms)')
ylabel('membrane potential (mV)')
title('voltage - 1c - 5ms delay')
grid on
saveas(f3, sprintf('1_c-5ms delay.png'));

applied_current = zeros(size(time_vector));
duration=5e-3/dt;
delay=10e-3/dt;
for j=1:10
    startpoint=((100e-3/dt)+delay*(j-1));
    endpoint=(startpoint+duration);
    applied_current(startpoint+1: endpoint) = 0.22e-9;
    
   

end
v = zeros(size(time_vector));
v(1) = leak_potential;
m = zeros(size(time_vector));
h = zeros(size(time_vector));
n = zeros(size(time_vector));
m(1) = 0;
h(1) = 0;
n(1) = 0;
for j = 1:length(time_vector)-1
    v(j+1) = v(j) + dt*( G_Leak*(leak_potential-v(j)) + G_Na*(m(j)^3)*h(j)*(E_Na-v(j))+G_K*(n(j)^4)*(E_K-v(j))+ applied_current(j))/C_m;
    if ( v(j) == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m=(1e5*(-v(j)-0.045))/(exp(100*(-v(j)-0.045))-1);
    end
    beta_m=4e3*exp((-v(j)-0.070)/0.018);
    alpha_h=70*exp(50*(-v(j)-0.070));
    beta_h=1e3/(1+exp(100*(-v(j)-0.040)));
    if ( v(j) == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else;                   % potassium activation rate
        alpha_n=(1e4*(-v(j)-0.060))/(exp(100*(-v(j)-0.060))-1);
    end
    beta_n=125*exp((-v(j)-0.070)/0.08);
    m(j+1)=m(j)+dt*(alpha_m*(1-m(j))-beta_m*m(j) );
    h(j+1)=h(j)+dt*(alpha_h*(1-h(j))-beta_h*h(j));
    n(j+1)=n(j)+dt*(alpha_n*(1-n(j))-beta_n*n(j));
end

f4=figure;
subplot(2,1,1);
plot(time_vector*1e3, applied_current*1e9)
xlabel('Time (ms)')
ylabel('Applied Current (nA)')
title('Applied Current with Step - 1c - 10ms delay')
grid on
subplot(2,1,2);
plot(time_vector*1e3, v*1e3)
xlabel('Time (ms)')
ylabel('membrane potential (mV)')
title('voltage - 1c - 10ms delay')
grid on
saveas(f4, sprintf('1_c-10ms delay.png'));

applied_current = zeros(size(time_vector));
duration=5e-3/dt;
delay=15e-3/dt;
for j=1:10
    startpoint=((100e-3/dt)+delay*(j-1));
    endpoint=(startpoint+duration);
    applied_current(startpoint+1: endpoint) = 0.22e-9;
    
   

end
v = zeros(size(time_vector));
v(1) = leak_potential;
m = zeros(size(time_vector));
h = zeros(size(time_vector));
n = zeros(size(time_vector));
m(1) = 0;
h(1) = 0;
n(1) = 0;
for j = 1:length(time_vector)-1
    v(j+1) = v(j) + dt*( G_Leak*(leak_potential-v(j)) + G_Na*(m(j)^3)*h(j)*(E_Na-v(j))+G_K*(n(j)^4)*(E_K-v(j))+ applied_current(j))/C_m;
    if ( v(j) == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m=(1e5*(-v(j)-0.045))/(exp(100*(-v(j)-0.045))-1);
    end
    beta_m=4e3*exp((-v(j)-0.070)/0.018);
    alpha_h=70*exp(50*(-v(j)-0.070));
    beta_h=1e3/(1+exp(100*(-v(j)-0.040)));
    if ( v(j) == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else;                   % potassium activation rate
        alpha_n=(1e4*(-v(j)-0.060))/(exp(100*(-v(j)-0.060))-1);
    end
    beta_n=125*exp((-v(j)-0.070)/0.08);
    m(j+1)=m(j)+dt*(alpha_m*(1-m(j))-beta_m*m(j) );
    h(j+1)=h(j)+dt*(alpha_h*(1-h(j))-beta_h*h(j));
    n(j+1)=n(j)+dt*(alpha_n*(1-n(j))-beta_n*n(j));
end

f5=figure;
subplot(2,1,1);
plot(time_vector*1e3, applied_current*1e9)
xlabel('Time (ms)')
ylabel('Applied Current (nA)')
title('Applied Current with Step - 1c - 15ms delay')
grid on
subplot(2,1,2);
plot(time_vector*1e3, v*1e3)
xlabel('Time (ms)')
ylabel('membrane potential (mV)')
title('voltage - 1c - 15ms delay')
grid on
saveas(f5, sprintf('1_c-15ms delay.png'));

applied_current = zeros(size(time_vector));
duration=5e-3/dt;
delay=25e-3/dt;
for j=1:10
    startpoint=((100e-3/dt)+delay*(j-1));
    endpoint=(startpoint+duration);
    applied_current(startpoint+1: endpoint) = 0.22e-9;
    
   

end
v = zeros(size(time_vector));
v(1) = leak_potential;
m = zeros(size(time_vector));
h = zeros(size(time_vector));
n = zeros(size(time_vector));
m(1) = 0;
h(1) = 0;
n(1) = 0;
for j = 1:length(time_vector)-1
    v(j+1) = v(j) + dt*( G_Leak*(leak_potential-v(j)) + G_Na*(m(j)^3)*h(j)*(E_Na-v(j))+G_K*(n(j)^4)*(E_K-v(j))+ applied_current(j))/C_m;
    if ( v(j) == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m=(1e5*(-v(j)-0.045))/(exp(100*(-v(j)-0.045))-1);
    end
    beta_m=4e3*exp((-v(j)-0.070)/0.018);
    alpha_h=70*exp(50*(-v(j)-0.070));
    beta_h=1e3/(1+exp(100*(-v(j)-0.040)));
    if ( v(j) == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else;                   % potassium activation rate
        alpha_n=(1e4*(-v(j)-0.060))/(exp(100*(-v(j)-0.060))-1);
    end
    beta_n=125*exp((-v(j)-0.070)/0.08);
    m(j+1)=m(j)+dt*(alpha_m*(1-m(j))-beta_m*m(j) );
    h(j+1)=h(j)+dt*(alpha_h*(1-h(j))-beta_h*h(j));
    n(j+1)=n(j)+dt*(alpha_n*(1-n(j))-beta_n*n(j));
end

f51=figure;
subplot(2,1,1);
plot(time_vector*1e3, applied_current*1e9)
xlabel('Time (ms)')
ylabel('Applied Current (nA)')
title('Applied Current with Step - 1c - 15ms delay')
grid on
subplot(2,1,2);
plot(time_vector*1e3, v*1e3)
xlabel('Time (ms)')
ylabel('membrane potential (mV)')
title('voltage - 1c - 15ms delay')
grid on
saveas(f51, sprintf('1_c-25ms delay.png'));
%Below a certain threshold: The delay might be too short for the sodium channels to recover from inactivation, preventing spike generation.
%Optimal range: There might be a "sweet spot" where the delay allows enough sodium channel recovery but keeps the membrane potential sufficiently depolarized for a strong response.
%Beyond a certain threshold: Very long delays might allow the membrane potential to return too close to the resting potential, making it less likely to reach the threshold even with high-amplitude, long-duration pulses
%Increased delay means less sodium channels are ready at the pulse onset, as they experience inactivation during the resting period. This initially reduces the inward sodium current, making it harder to reach the threshold.
%However, a longer delay also allows more time for the potassium channels to close, increasing the overall membrane resistance and promoting further depolarization.
%Depending on the interplay of these opposing effects, increasing the delay can sometimes lead to larger spikes due to greater membrane resistance if the sodium channels recover enough by the time the pulse arrives.
%Very long delays typically result in no spikes as the membrane potential returns closer to resting potential, requiring a stronger or longer pulse to overcome the leak current.
%Time-dependent inactivation: Sodium channels have an intrinsic inactivation mechanism. Once activated by a depolarization, they deactivate and enter an inactivated state over time. The recovery from this inactivation also takes time, following a specific rate constant.
%Leak current: During the resting period between pulses, the leak current counteracts depolarization and gradually pulls the membrane potential back towards its resting level.
%Interplay of recovery and leak: While a longer delay (like 15 ms) offers more time for recovery, it also allows the leak current more time to repolarize the membrane. This can counteract the sodium channel recovery, making it harder for the membrane potential to reach the threshold and trigger a spike.
%Why 10 ms might be an "optimal" delay in your case:

%With a 10 ms delay, the sodium channels might partially recover, but not completely. This can be enough to generate a significant inward current upon pulse arrival, pushing the membrane potential closer to the threshold.
%Additionally, the 10 ms delay might not give the leak current enough time to fully repolarize the membrane, leaving it still sufficiently depolarized for the recovered sodium channels to contribute a strong response.
%Why 15 ms might not trigger a spike:

%As the delay increases to 15 ms, the leak current has more time to act, potentially pulling the membrane potential further away from the threshold.
%Even though some sodium channels might recover from inactivation, the initial membrane potential could be too low to overcome the leak current and the remaining inactivated channels, preventing a strong enough depolarization to reach the threshold.

%% part D
applied_current = 0.6e-9*ones(size(time_vector));
v = zeros(size(time_vector));
v(1) = -0.065;
m = zeros(size(time_vector));
h = zeros(size(time_vector));
n = zeros(size(time_vector));
m(1) = 0.05;
h(1) = 0.5;
n(1) = 0.35;
duration=5e-3/dt;
delay=0;
for j=1:10
    startpoint=((100e-3/dt)+delay+duration*(j-1));
    endpoint=(startpoint+duration);
    applied_current(startpoint: endpoint) =0;
    delay=j*20e-3/dt;
    
end
for j = 1:length(time_vector)-1
    v(j+1) = v(j) + dt*( G_Leak*(leak_potential-v(j)) + G_Na*(m(j)^3)*h(j)*(E_Na-v(j))+G_K*(n(j)^4)*(E_K-v(j))+ applied_current(j))/C_m;
    if ( v(j) == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m=(1e5*(-v(j)-0.045))/(exp(100*(-v(j)-0.045))-1);
    end
    beta_m=4e3*exp((-v(j)-0.070)/0.018);
    alpha_h=70*exp(50*(-v(j)-0.070));
    beta_h=1e3/(1+exp(100*(-v(j)-0.040)));
    if ( v(j) == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else;                   % potassium activation rate
        alpha_n=(1e4*(-v(j)-0.060))/(exp(100*(-v(j)-0.060))-1);
    end
    beta_n=125*exp((-v(j)-0.070)/0.08);
    m(j+1)=m(j)+dt*(alpha_m*(1-m(j))-beta_m*m(j) );
    h(j+1)=h(j)+dt*(alpha_h*(1-h(j))-beta_h*h(j));
    n(j+1)=n(j)+dt*(alpha_n*(1-n(j))-beta_n*n(j));
end
f6=figure;
subplot(2,1,1);
plot(time_vector*1e3, applied_current*1e9)
xlabel('Time (ms)')
ylabel('Applied Current (nA)')
title('Applied Current with Step - 1d ')
grid on
subplot(2,1,2);
plot(time_vector*1e3, v*1e3)
xlabel('Time (ms)')
ylabel('membrane potential (mV)')
title('voltage - 1d')
grid on
saveas(f6, sprintf('1_d.png'));

%I see that

%% part e
applied_current = 0.65e-9*ones(size(time_vector));
v = zeros(size(time_vector));
v(1) = -0.065;
m = zeros(size(time_vector));
h = zeros(size(time_vector));
n = zeros(size(time_vector));
m(1) = 0.05;
h(1) = 0.5;
n(1) = 0.35;
duration=5e-3/dt;
delay=0;
for j=1
    startpoint=((100e-3/dt)+delay+duration*(j-1));
    endpoint=(startpoint+duration);
    applied_current(startpoint: endpoint) =1e-9;
    delay=j*20e-3/dt;
   
end
for j = 1:length(time_vector)-1
    v(j+1) = v(j) + dt*( G_Leak*(leak_potential-v(j)) + G_Na*(m(j)^3)*h(j)*(E_Na-v(j))+G_K*(n(j)^4)*(E_K-v(j))+ applied_current(j))/C_m;
    if ( v(j) == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m=(1e5*(-v(j)-0.045))/(exp(100*(-v(j)-0.045))-1);
    end
    G_Na*(m(j)^3)*h(j)*(E_Na-v(j))
    beta_m=4e3*exp((-v(j)-0.070)/0.018);
    alpha_h=70*exp(50*(-v(j)-0.070));
    beta_h=1e3/(1+exp(100*(-v(j)-0.040)));
    if ( v(j) == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else;                   % potassium activation rate
        alpha_n=(1e4*(-v(j)-0.060))/(exp(100*(-v(j)-0.060))-1);
    end
    G_K*(n(j)^4)*(E_K-v(j))
    beta_n=125*exp((-v(j)-0.070)/0.08);
    m(j+1)=m(j)+dt*(alpha_m*(1-m(j))-beta_m*m(j) );
    h(j+1)=h(j)+dt*(alpha_h*(1-h(j))-beta_h*h(j));
    n(j+1)=n(j)+dt*(alpha_n*(1-n(j))-beta_n*n(j));
end

f7=figure;
subplot(2,1,1);
plot(time_vector*1e3, applied_current*1e9)
xlabel('Time (ms)')
ylabel('Applied Current (nA)')
title('Applied Current with Step - 1e ')
grid on
subplot(2,1,2);
plot(time_vector*1e3, v*1e3)
xlabel('Time (ms)')
ylabel('membrane potential (mV)')
title('voltage - 1e')
grid on
saveas(f7, sprintf('1_e.png'));
%% part f
applied_current = 0.7e-9*ones(size(time_vector));
v = zeros(size(time_vector));
v(1) = -0.065;
m = zeros(size(time_vector));
h = zeros(size(time_vector));
n = zeros(size(time_vector));
m(1) = 0;
h(1) =0;
n(1) = 0;
duration=5e-3/dt;
delay=0;
for j=1
    startpoint=((100e-3/dt)+delay+duration*(j-1));
    endpoint=(startpoint+duration);
    applied_current(startpoint: endpoint) =1e-9;
    delay=j*20e-3/dt;
   
end
for j = 1:length(time_vector)-1
    v(j+1) = v(j) + dt*( G_Leak*(leak_potential-v(j)) + G_Na*(m(j)^3)*h(j)*(E_Na-v(j))+G_K*(n(j)^4)*(E_K-v(j))+ applied_current(j))/C_m;
    if ( v(j) == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m=(1e5*(-v(j)-0.045))/(exp(100*(-v(j)-0.045))-1);
    end
    beta_m=4e3*exp((-v(j)-0.070)/0.018);
    alpha_h=70*exp(50*(-v(j)-0.070));
    beta_h=1e3/(1+exp(100*(-v(j)-0.040)));
    if ( v(j) == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else;                   % potassium activation rate
        alpha_n=(1e4*(-v(j)-0.060))/(exp(100*(-v(j)-0.060))-1);
    end
    beta_n=125*exp((-v(j)-0.070)/0.08);
    m(j+1)=m(j)+dt*(alpha_m*(1-m(j))-beta_m*m(j) );
    h(j+1)=h(j)+dt*(alpha_h*(1-h(j))-beta_h*h(j));
    n(j+1)=n(j)+dt*(alpha_n*(1-n(j))-beta_n*n(j));
end
f8=figure;
subplot(2,1,1);
plot(time_vector*1e3, applied_current*1e9)
xlabel('Time (ms)')
ylabel('Applied Current (nA)')
title('Applied Current with Step - 1f ')
grid on
subplot(2,1,2);
plot(time_vector*1e3, v*1e3)
xlabel('Time (ms)')
ylabel('membrane potential (mV)')
title('voltage - 1f')
grid on
saveas(f8, sprintf('1_f.png'));
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');