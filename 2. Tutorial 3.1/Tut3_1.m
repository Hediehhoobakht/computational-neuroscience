close all
clear 
clc
 
leak_potential = -60e-3;               % Leak potential (V)
V_threshold = -50e-3;         % Threshold potential (V)
reset_potential = -80e-3;          % Reset potential (V)
delta_th = 2e-3;             % Threshold shift factor (V)
G_Leak = 8e-9;                % Leak conductance (S)
R=1/8e9;                    % resistance (ohm)
C = 100e-12;                % Capacitance (F)
a = 10e-9;                  % adaptation recovery (S)
b = 0.5e-9;                % adaptation strength (A)
tau_SRA = 50e-3;             % Adaptation time constant (s)

V_max = 50e-3;              % level of voltage to detect a spike
tau_m = R*C;

current_val = 1e-9 * (rand(1, 40000) -0.5);  % Scale to nA and shift to range [-0.5, 0.5]

total_time = 200000e-3;
dt=0.02e-3;
time_vector = 0:dt:total_time;

applied_current = zeros(size(time_vector));
for j=1:(length(current_val))
for i=1+250*(j-1):250*j
    applied_current(i)=current_val(j);
end
end
v = zeros(size(time_vector));       % initialize voltage
v(1) = leak_potential;
I_sra = zeros(size(time_vector));       % initialize adaptation variable
spikes = zeros(size(time_vector));

for j = 1:length(time_vector)-1     

    if ( v(j) > V_max )          % if there is a spike
        v(j) = reset_potential;         % reset the voltage
        I_sra(j) = I_sra(j) + b;        % increase the adaptation variable by b
        spikes(j) = 1;          % Record the spike as a "1" in the time-bin
    end
    
    % next line integrates the voltage over time, first part is like LIF
    % second part is an exponential spiking term
    % third part includes adaptation
    v(j+1) = v(j) + dt*( G_Leak*(leak_potential-v(j) + delta_th*exp((v(j)-V_threshold)/delta_th) ) ...
       - I_sra(j) + applied_current(j))/C;

   % next line decys the adaptation toward a steady state in between spikes
    I_sra(j+1) = I_sra(j) + dt*( a*(v(j)-leak_potential) - I_sra(j) )/tau_SRA;
    
end
f1 = figure;
figure(f1);
plot(time_vector, v)
grid on
title('Membrane Potential vs. Time');
xlabel('Time');
ylabel('Membrane Potential');   
saveas(f1, sprintf('membrane_potential.png'));

% Now downsample the stimulus and response to 1ms bins using the online
% function expandbin
new_dt = 0.001;                      % new time-bin of 1ms
spikes = expandbin(spikes,dt,new_dt);
spikes(find(spikes)) = 1;   % replace fractions with a "1" for a spike
applied_current = expandbin(applied_current,dt,new_dt);  % The function expandbin does the downsampling
new_time_vector = 0:new_dt:total_time;

f2 = figure;
plot(new_time_vector(1:200000), applied_current(1:200000))
grid on
title('applied_current vs. Time');
xlabel('Time');
ylabel('I_{app}');
saveas(f2, sprintf('applied_current.png'));
[sta, tcorr] = STA(applied_current, spikes, new_dt);
f3 = figure;
plot(tcorr,sta)
xlabel('Spike lag (ms)')
ylabel('Stimulus strength')
grid on
title('Stimulus_strength vs. Spike lag');
saveas(f3, sprintf('Stimulus_strength.png'));

%%%%%%%%%%%%
%B
Nsteps = 40000;             % Number of time-steps of distinct applied currents
Nspatial = 40;              % Number of spatially distinct values per time-step
s = (rand(Nspatial, Nsteps) - 0.5)* 1e-9;

x_max = 40; % Assuming x_max value
x_0 = 20.5;    % Assuming x_0 value
I_app = zeros(1, Nsteps);
for t = 1:length(I_app)-1
    sum_val = 0;
    for x = 1:x_max
        sum_val = sum_val + cos(4*pi*(x - x_0)/x_max) .* exp(-16*((x - x_0)/x_max).^2) .* s(x, t); % Assuming s(x, t) is a function defined elsewhere
    end
    I_app(t) = sum_val;
end
x = 1:x_max;
w=cos(4*pi*(x - x_0)/x_max) .* exp(-16*((x - x_0)/x_max).^2);
% Plot the weight vector
f4=figure;
plot(x, w);
xlabel('Spatial coordinate x');
ylabel('Weight w(x)');
title('Input weight vector w(x)');
grid on
saveas(f4, sprintf('Spatial coordinate vs Input weight vector_w.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
a = 40e-9;  
b = 1e-9;
total_time = 200000e-3;
dt = 0.02e-3;
time_vector = 0:dt:total_time;

I_app_B = zeros(Nspatial,length(time_vector));
for j=1:(length(I_app))
for i=1+250*(j-1):250*j
    I_app_B(:,i)=I_app(j);
end
end
% Initialize variables
v = zeros(size(time_vector));
v(1) = leak_potential;
I_sra = zeros(size(time_vector));
spikes_B = zeros(size(time_vector));

% Simulation loop
for j = 1:length(time_vector)-1

    if (v(j) > V_max)
        v(j) = reset_potential;
        I_sra(j) = I_sra(j) + b;
        spikes_B(j) = 1;
    end

    v(j+1) = v(j) + dt*(G_Leak*(leak_potential-v(j) + delta_th*exp((v(j)-V_threshold)/delta_th)) ...
            - I_sra(j) + I_app_B(j))/C;

    I_sra(j+1) = I_sra(j) + dt*(a*(v(j)-leak_potential) - I_sra(j))/tau_SRA;

end

% Plot membrane potential
f5=figure;
plot(time_vector, v)
grid on
title('Membrane Potential vs. Time');
xlabel('Time');
ylabel('Membrane Potential');
saveas(f5, sprintf('Membrane Potential vs. Time part_B.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_dt = 0.001;                      % new time-bin of 1ms
spikes_B = expandbin(spikes_B,dt,new_dt);
spikes_B(find(spikes_B)) = 1;

%%%%%%%%%%%%%%%
steplength = 0.005;
new_Nt = length(spikes_B);                 % New length of spike vector
newIapp = zeros(Nspatial,new_Nt);        % Define a new input vector
new_nstep_length = round(steplength/new_dt); % New number of 1ms bins per step

for step = 1:Nsteps;                    % Loop through steps of constant current
    istart = (step-1)*new_nstep_length+1; % first time point with 1ms steps
    istop = step*new_nstep_length;        % last time point with 1ms steps
    
    % generate the input vector as before, but with bins of 1ms, not of dt
    new_Iapp(:,istart:istop) = s(:,step)*ones(1,new_nstep_length);
end
%%%%%%%%%%


[sta, tcorr] = STA_spatial(new_Iapp, spikes_B, new_dt,new_nstep_length); 

f6=figure()
imagesc(fliplr(sta));       % reverses time-axis to plot STA
colormap(gray)              % grayscale
set(gca,'XTick',[1, 26, 51, 76, 101])
set(gca,'XTickLabel',{'-25' '0', '25', '50', '75'})
xlabel('Spike lag (ms)')
ylabel('Stimulus coordinate')
set(gca,'YTick',[ ])
saveas(f6, sprintf('imagsec.png'));
f7=figure();
plot(sta(12, :));title('STA row12');saveas(f7, sprintf('row12.png'));
f8=figure();
plot(sta(20, :));title('STA row20');saveas(f8, sprintf('row20.png'));
f9=figure();
plot(sta(28, :));title('STA row28');saveas(f9, sprintf('row28.png'));
f10=figure();
plot(sta(:, 25));title('STA col25');saveas(f10, sprintf('col25.png'));
f11=figure();
plot(sta(:, 50));title('STA col50');saveas(f11, sprintf('col50.png'));
f12=figure();
plot(sta(:, 75));title('STA col75');saveas(f12, sprintf('col75.png'));
