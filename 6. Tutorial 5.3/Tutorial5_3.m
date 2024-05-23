close all
clear
clc
%define parameters
dt=1e-3; %I consider the time interval of 1 ms
t=0:dt:6-dt;
I_0 = 2e-9;
I_app1 = I_0*ones(size(t));
I_app2 = I_0*ones(size(t));
V1 = zeros(size(t));
V2 = zeros(size(t));
s1 = zeros(size(t));
s2 = zeros(size(t));
C_m = 1e-9;
R_m = 10e6;
tau_m=C_m*R_m;
E_L = -70e-3;
Vth = -54e-3;
Vreset = -80e-3;
Erev_12 = -70e-3;
Erev_21 = -70e-3;
G_12 = 1e-6;
G_21 = 1e-6;
tau_syn = 10e-3;

p_R = 1;
D1(1) = 1;
D2(1) = 1;
s1(1) = 0;
s2(1) = 0.;
tau_D=0.2;


V1(1) = E_L;
V2(1) = E_L;
%% A2
sigma_I=0;
for n = 2:length(t)
    if n<=0.1/dt+1;
        I_app1(n-1) = I_app1(n-1) + 3e-9;

    end
    if 30*0.1/dt<=n && n<=0.1/dt+30*0.1/dt-1;
        I_app2(n-1) =I_app2(n-1) + 3e-9;

    end
   V1(n)= V1(n-1) + ((E_L - V1(n-1))/tau_m +G_21*s2(n-1)*(Erev_21-V1(n-1))/C_m + I_app1(n-1)/C_m +randn*sigma_I/C_m)*dt;
   s1(n) = s1(n-1) - dt*s1(n-1)/tau_syn;
   D1(n) = D1(n-1)+dt*(1-D1(n-1))/tau_D;
   V2(n)= V2(n-1) + ((E_L - V2(n-1))/tau_m +G_12*s1(n-1)*(Erev_12-V2(n-1))/C_m + I_app2(n-1)/C_m +randn*sigma_I/C_m)*dt;
   s2(n) = s2(n-1) - dt*s2(n-1)/tau_syn;
   D2(n) = D2(n-1)+dt*(1-D2(n-1))/tau_D;
   if V1(n) > Vth
      V1(n) = Vreset;
      s1(n) = s1(n) + p_R*D1(n)*(1-s1(n));
      D1(n)=D1(n)*(1-p_R);
   end
   if V2(n) > Vth
      V2(n) = Vreset;
      s2(n) = s2(n) + p_R*D2(n)*(1-s2(n));
      D2(n)=D2(n)*(1-p_R);
   end
end
% Plot results
f1=figure(1);
subplot(3,1,1);
plot(t*1e3, V1, 'b', t*1e3, V2, 'r');
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
legend('Neuron 1', 'Neuron 2');
title('Membrane potential of the two neurons');

subplot(3,1,2);
plot(t*1e3, s1, 'b', t*1e3, s2, 'r');
xlabel('Time (ms)');
ylabel('Synaptic gating variable');
legend('s1','s2');
title('Gating variables of the two neurons');

subplot(3,1,3);
plot(t*1e3, I_app1*1e9, 'b', t*1e3, I_app2*1e9, 'r');
xlabel('Time (ms)');
ylabel('I-app (nA)');
legend('I-app1','I-app2');
title('applied current of the two neurons');
saveas(f1, sprintf('A2.png'));

%% A3

I_app1 = I_0*ones(size(t));
I_app2 = I_0*ones(size(t));
V1 = zeros(size(t));
V2 = zeros(size(t));
s1 = zeros(size(t));
s2 = zeros(size(t));
V1(1) = E_L;
V2(1) = E_L;
D1(1) = 1;
D2(1) = 1;
s1(1) = 0;
s2(1) = 0.;
sigma_I=50e-12;
for n = 2:length(t)
   V1(n)= V1(n-1) + ((E_L - V1(n-1))/tau_m +G_21*s2(n-1)*(Erev_21-V1(n-1))/C_m + I_app1(n-1)/C_m +randn*sigma_I/C_m)*dt;
   s1(n) = s1(n-1) - dt*s1(n-1)/tau_syn;
   D1(n) = D1(n-1)+dt*(1-D1(n-1))/tau_D;
   V2(n)= V2(n-1) + ((E_L - V2(n-1))/tau_m +G_12*s1(n-1)*(Erev_12-V2(n-1))/C_m + I_app2(n-1)/C_m +randn*sigma_I/C_m)*dt;
   s2(n) = s2(n-1) - dt*s2(n-1)/tau_syn;
   D2(n) = D2(n-1)+dt*(1-D2(n-1))/tau_D;
   if V1(n) > Vth
      V1(n) = Vreset;
      s1(n) = s1(n) + p_R*D1(n)*(1-s1(n));
      D1(n)=D1(n)*(1-p_R);
   end
   if V2(n) > Vth
      V2(n) = Vreset;
      s2(n) = s2(n) + p_R*D2(n)*(1-s2(n));
      D2(n)=D2(n)*(1-p_R);
   end
end
% Plot results
f2=figure(2);
subplot(3,1,1);
plot(t*1e3, V1, 'b', t*1e3, V2, 'r');
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
legend('Neuron 1', 'Neuron 2');
title('Membrane potential of the two neurons');

subplot(3,1,2);
plot(t*1e3, s1, 'b', t*1e3, s2, 'r');
xlabel('Time (ms)');
ylabel('Synaptic gating variable');
legend('s1','s2');
title('Gating variables of the two neurons');

subplot(3,1,3);
plot(t*1e3, I_app1*1e9, 'b', t*1e3, I_app2*1e9, 'r');
xlabel('Time (ms)');
ylabel('I-app (nA)');
legend('I-app1','I-app2');
title('applied current of the two neurons');
saveas(f2, sprintf('A3.png'));

%% A4
dt=1e-3; %I consider the time interval of 1 ms
t=0:dt:6-dt;
I_app1 = I_0*ones(size(t));
I_app2 = I_0*ones(size(t));
V1 = zeros(size(t));
V2 = zeros(size(t));
s1 = zeros(size(t));
s2 = zeros(size(t));
V1(1) = E_L;
V2(1) = E_L;
state=1;
switch_times=[];
num_switches=0;
sigma_I=50e-12;
D1(1) = 1;
D2(1) = 1;
s1(1) = 0;
s2(1) = 0.;
for n = 2:length(t)
   V1(n)= V1(n-1) + ((E_L - V1(n-1))/tau_m +G_21*s2(n-1)*(Erev_21-V1(n-1))/C_m + I_app1(n-1)/C_m +randn*sigma_I/C_m)*dt;
   s1(n) = s1(n-1) - dt*s1(n-1)/tau_syn;
   D1(n) = D1(n-1)+dt*(1-D1(n-1))/tau_D;
   V2(n)= V2(n-1) + ((E_L - V2(n-1))/tau_m +G_12*s1(n-1)*(Erev_12-V2(n-1))/C_m + I_app2(n-1)/C_m +randn*sigma_I/C_m)*dt;
   s2(n) = s2(n-1) - dt*s2(n-1)/tau_syn;
   D2(n) = D2(n-1)+dt*(1-D2(n-1))/tau_D;
   if V1(n) > Vth && state==2
      V1(n) = Vreset;
      s1(n) = s1(n) + p_R*D1(n)*(1-s1(n));
      D1(n)=D1(n)*(1-p_R);
      state=1;
      num_switches=num_switches+1;
      switch_times=[switch_times,t(n)];
   end
   if V2(n) > Vth && state==1;
      V2(n) = Vreset;
      s2(n) = s2(n) + p_R*D2(n)*(1-s2(n));
      D2(n)=D2(n)*(1-p_R);
      state=2;
      num_switches=num_switches+1;
      switch_times=[switch_times,t(n)];
   end
end
num_switches
% Plot results
f3=figure(3);
subplot(3,1,1);
plot(t*1e3, V1, 'b', t*1e3, V2, 'r');
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
legend('Neuron 1', 'Neuron 2');
title('Membrane potential of the two neurons');

subplot(3,1,2);
plot(t*1e3, s1, 'b', t*1e3, s2, 'r');
xlabel('Time (ms)');
ylabel('Synaptic gating variable');
legend('s1','s2');
title('Gating variables of the two neurons');

subplot(3,1,3);
plot(t*1e3, I_app1*1e9, 'b', t*1e3, I_app2*1e9, 'r');
xlabel('Time (ms)');
ylabel('I-app (nA)');
legend('I-app1','I-app2');
title('applied current of the two neurons');
saveas(f3, sprintf('A4.png'));
durations = diff(switch_times);
f4=figure;
subplot(2,1,1);
histogram(durations(1:2:end), 'FaceColor', 'b');
title('State 1 durations');
xlabel('Duration (ms)');
ylabel('Frequency');

subplot(2,1,2);
histogram(durations(2:2:end), 'FaceColor', 'r');
title('State 2 durations');
xlabel('Duration (ms)');

saveas(f4, 'A4-hist.png');
%% b1
dt=1e-3; %I consider the time interval of 1 ms
t=0:dt:6-dt;
p_R = 0.2;
V1 = zeros(size(t));
V2 = zeros(size(t));
s1 = zeros(size(t));
s2 = zeros(size(t));
D1(1) = 1;
D2(1) = 1;
s1(1) = 0;
s2(1) = 0.;
tau_D=0.2;
sigma_I=0;
I_app1 = I_0*ones(size(t));
I_app2 = I_0*ones(size(t));

V1(1) = E_L;
V2(1) = E_L;
for n = 2:length(t)
    if n<=0.1/dt+1;
        I_app1(n-1) = I_app1(n-1) + 3e-9;

    end
    if 30*0.1/dt<=n && n<=0.1/dt+30*0.1/dt-1;
        I_app2(n-1) =I_app2(n-1) + 3e-9;

    end
   V1(n)= V1(n-1) + ((E_L - V1(n-1))/tau_m +G_21*s2(n-1)*(Erev_21-V1(n-1))/C_m + I_app1(n-1)/C_m +randn*sigma_I/C_m)*dt;
   s1(n) = s1(n-1) - dt*s1(n-1)/tau_syn;
   D1(n) = D1(n-1)+dt*(1-D1(n-1))/tau_D;
   V2(n)= V2(n-1) + ((E_L - V2(n-1))/tau_m +G_12*s1(n-1)*(Erev_12-V2(n-1))/C_m + I_app2(n-1)/C_m +randn*sigma_I/C_m)*dt;
   s2(n) = s2(n-1) - dt*s2(n-1)/tau_syn;
   D2(n) = D2(n-1)+dt*(1-D2(n-1))/tau_D;
   if V1(n) > Vth
      V1(n) = Vreset;
      s1(n) = s1(n) + p_R*D1(n)*(1-s1(n));
      D1(n)=D1(n)*(1-p_R);
   end
   if V2(n) > Vth
      V2(n) = Vreset;
      s2(n) = s2(n) + p_R*D2(n)*(1-s2(n));
      D2(n)=D2(n)*(1-p_R);
   end
end
% Plot results
f5=figure(5);
subplot(3,1,1);
plot(t*1e3, V1, 'b', t*1e3, V2, 'r');
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
legend('Neuron 1', 'Neuron 2');
title('Membrane potential of the two neurons');

subplot(3,1,2);
plot(t*1e3, s1, 'b', t*1e3, s2, 'r');
xlabel('Time (ms)');
ylabel('Synaptic gating variable');
legend('s1','s2');
title('Gating variables of the two neurons');

subplot(3,1,3);
plot(t*1e3, I_app1*1e9, 'b', t*1e3, I_app2*1e9, 'r');
xlabel('Time (ms)');
ylabel('I-app (nA)');
legend('I-app1','I-app2');
title('applied current of the two neurons');
saveas(f5, sprintf('b1.png'));
%% b2
I_app1 = I_0*ones(size(t));
I_app2 = I_0*ones(size(t));
V1 = zeros(size(t));
V2 = zeros(size(t));
s1 = zeros(size(t));
s2 = zeros(size(t));
D1(1) = 1;
D2(1) = 1;
s1(1) = 0;
s2(1) = 0.;

V1(1) = E_L;
V2(1) = E_L;
state=1;
switch_times=[];
num_switches=0;
sigma_I=5e-12;
for n = 2:length(t)
   V1(n)= V1(n-1) + ((E_L - V1(n-1))/tau_m +G_21*s2(n-1)*(Erev_21-V1(n-1))/C_m + I_app1(n-1)/C_m +randn*sigma_I/C_m)*dt;
   s1(n) = s1(n-1) - dt*s1(n-1)/tau_syn;
   D1(n) = D1(n-1)+dt*(1-D1(n-1))/tau_D;
   V2(n)= V2(n-1) + ((E_L - V2(n-1))/tau_m +G_12*s1(n-1)*(Erev_12-V2(n-1))/C_m + I_app2(n-1)/C_m +randn*sigma_I/C_m)*dt;
   s2(n) = s2(n-1) - dt*s2(n-1)/tau_syn;
   D2(n) = D2(n-1)+dt*(1-D2(n-1))/tau_D;
   if V1(n) > Vth && state==2
      V1(n) = Vreset;
      s1(n) = s1(n) + p_R*D1(n)*(1-s1(n));
      D1(n)=D1(n)*(1-p_R);
      state=1;
      num_switches=num_switches+1;
      switch_times=[switch_times,t(n)];
   end
   if V2(n) > Vth && state==1;
      V2(n) = Vreset;
      s2(n) = s2(n) + p_R*D2(n)*(1-s2(n));
      D2(n)=D2(n)*(1-p_R);
      state=2;
      num_switches=num_switches+1;
      switch_times=[switch_times,t(n)];
   end
end
% Plot results
f6=figure(6);
subplot(3,1,1);
plot(t*1e3, V1, 'b', t*1e3, V2, 'r');
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
legend('Neuron 1', 'Neuron 2');
title('Membrane potential of the two neurons');

subplot(3,1,2);
plot(t*1e3, s1, 'b', t*1e3, s2, 'r');
xlabel('Time (ms)');
ylabel('Synaptic gating variable');
legend('s1','s2');
title('Gating variables of the two neurons');

subplot(3,1,3);
plot(t*1e3, I_app1*1e9, 'b', t*1e3, I_app2*1e9, 'r');
xlabel('Time (ms)');
ylabel('I-app (nA)');
legend('I-app1','I-app2');
title('applied current of the two neurons');
saveas(f6, sprintf('b2.png'));
durations = diff(switch_times);
f7=figure;
subplot(2,1,1);
histogram(durations(1:2:end), 'FaceColor', 'b');
title('State 1 durations');
xlabel('Duration (ms)');
ylabel('Frequency');

subplot(2,1,2);
histogram(durations(2:2:end), 'FaceColor', 'r');
title('State 2 durations');
xlabel('Duration (ms)');

saveas(f7, 'b2-hist.png');