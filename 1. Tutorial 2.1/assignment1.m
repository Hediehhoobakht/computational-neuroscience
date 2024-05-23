close all
clear all
clc
% Define parameters
C = 2e-9;  % Membrane capacitance (F)
R = 5*1e6;     % Membrane resistance (?)
leak_potential = -70e-3;  % Leak potential (V)
threshold = -50e-3;      % Spike threshold (V)
reset_potential = -65e-3;  % Reset potential (V)
dt = 0.1e-3;     % Time step (s)
T = 2;      % Maximum time (s)

% Create time vector
t = 0:dt:T;

% I define the I_app and the loop after that in Iapp function
I_th = (threshold  - leak_potential)/R;
I_lower = I_th*0.99;
[I_app_lower,v_lower] = Iapp(I_lower, t,leak_potential,C,R,dt,threshold,reset_potential);
I_higher = I_th*1.01;
[I_app_higher,v_higher] = Iapp(I_higher, t,leak_potential,C,R,dt,threshold,reset_potential);

img=figure;
plot(t(1:2001), v_lower(1:2001), 'b', 'LineWidth', 2);
hold on;
plot(t(1:2001), v_higher(1:2001), 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Membrane Potential (V)');
title('Membrane Potential with Different Applied Currents');
legend('slightly lower than','slightly higher than','Location','SouthEast');
grid on;
saveas(img, sprintf('Tutorial_2_1_question_1b.png'));
%Lower Current: I_low doesn't produce spikes, the calculation and simulation confirm the minimum current threshold for spiking.
%Higher Current:I_high produces spikes, it verifies that currents above the threshold trigger spiking behavior.
Tau = C*R;
i_0 = [I_th*1.05 I_th*1.1 I_th*1.15 I_th*1.2 I_th*1.25 I_th*1.30 I_th*1.25 I_th*1.35 I_th*1.4 I_th*1.45];
frates = zeros(size(i_0));
k=1;
for I_0 = i_0
    
    [I_app_f,v_f]  = Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential);
   
    spike_count = 0;
    
    for n = 1:length(v_f)
        if v_f(n) == reset_potential;
            spike_count = spike_count + 1;
        end
    end
    rate = spike_count/2;
    frates(k) = rate;
    k=k+1;
end
img1=figure
scatter(i_0, frates, 'LineWidth', 2);
xlabel('I_app');
ylabel('Firing rate(Hz)');
title('Firing rate based on injected current 1c');
legend('1c','Location','SouthEast');
grid on;
saveas(img1, sprintf('Tutorial_2_1_question_1c.png'));

f=zeros(size(i_0));
j=1;
for I_0 = i_0
    if ((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - reset_potential) > 0) & ...
       ((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - threshold) > 0)
        
        firing_rate = 1 / (Tau * log((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - reset_potential) / ...
                       (leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - threshold)));
    
    end
    
    f(j) = firing_rate;
    j = j + 1;
end

img2=figure
scatter(i_0, frates, 'LineWidth', 2);
hold on
scatter(i_0, f, 'LineWidth', 2);
xlabel('I_app');
ylabel('Firing rate(Hz)');
title('Firing rate based on injected current');
legend('1c','1d','Location','SouthEast');
grid on;
saveas(img2, sprintf('Tutorial_2_1_question_1d.png'));
sigma=[0.1,0.2]

I_lower = I_th*0.99;
[I_app_lower_n,v_lower_n] = n_Iapp(I_lower, t,leak_potential,C,R,dt,threshold,reset_potential,sigma(1));
I_higher = I_th*1.01;
[I_app_higher_n,v_higher_n] = n_Iapp(I_higher, t,leak_potential,C,R,dt,threshold,reset_potential,sigma(1));

img3=figure;
plot(t(1:2001), v_lower(1:2001), 'b', 'LineWidth', 2);
hold on;
plot(t(1:2001), v_higher(1:2001), 'r', 'LineWidth', 2);
hold on
plot(t(1:2001), v_lower_n(1:2001), 'g', 'LineWidth', 2);
hold on;
plot(t(1:2001), v_higher_n(1:2001), 'c', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Membrane Potential (V)');
title('Membrane Potential with Different Applied Currents in presense of noise');
legend('slightly lower than','slightly higher than','slightly lower than-with noise','slightly higher than-with noise','Location','SouthEast');
grid on;
saveas(img3, sprintf('Tutorial_2_1_question_2a_1.png'));


n_frates = zeros(length(i_0),length(sigma));
k=1;
l=1;
for sigma_I=sigma
    
for I_0 = i_0
    
    [I_app_f,v_f]  = n_Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential,sigma_I);
   
    spike_count = 0;
    
    for n = 1:length(v_f)
        if v_f(n) == reset_potential;
            spike_count = spike_count + 1;
        end
    end
    rate = spike_count/2;
    n_frates(k,l) = rate;
    k=k+1;
end
k=1;
l=l+1;
end

f_n=zeros(length(i_0),length(sigma));
j=1;
w=1;
for sigma_I=sigma
for I_0 = i_0
    if ((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - reset_potential) > 0) & ...
       ((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - threshold) > 0)
        
        firing_rate = 1 / (Tau * log((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - reset_potential) / ...
                       (leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - threshold)));
    
    end
    
    f_n(j,w) = firing_rate;
    j = j + 1;
end
j=1;
w=w+1
end

img4=figure
scatter(i_0, frates, 'SizeData', 52, 'LineWidth', 2);
hold on
scatter(i_0, f, 'SizeData', 52, 'LineWidth', 2);
hold on
scatter(i_0, n_frates(:,1), 'Marker', '*', 'SizeData', 52, 'LineWidth', 2);
hold on
scatter(i_0, f_n(:,1), 'Marker', '*', 'SizeData', 52, 'LineWidth', 2);
hold on
scatter(i_0, n_frates(:,2), 'Marker', 's', 'SizeData', 52, 'LineWidth', 2);
hold on
scatter(i_0, f_n(:,2), 'Marker', 's', 'SizeData',52, 'LineWidth', 2);
xlabel('I_app');
ylabel('Firing rate(Hz)');
title('Firing rate based on injected current');
legend('1c-sigma=0','1d-sigma=0','1c-sigma=0.1','1d-sigma=0.1','1c-sigma=0.2','1d-sigma=0.2','Location','SouthEast');
grid on;
saveas(img4, sprintf('Tutorial_2_1_question_2a.png'));
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.01e-3;     % Time step (s)
T = 2;      % Maximum time (s)

% Create time vector
t = 0:dt:T;

% I define the I_app and the loop after that in Iapp function
I_th = (threshold  - leak_potential)/R;
I_lower = I_th*0.99;
[I_app_lower,v_lower] = Iapp(I_lower, t,leak_potential,C,R,dt,threshold,reset_potential);
I_higher = I_th*1.01;
[I_app_higher,v_higher] = Iapp(I_higher, t,leak_potential,C,R,dt,threshold,reset_potential);

img5=figure;
plot(t(1:2001), v_lower(1:2001), 'b', 'LineWidth', 2);
hold on;
plot(t(1:2001), v_higher(1:2001), 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Membrane Potential (V)');
title('Membrane Potential with Different Applied Currents');
legend('slightly lower than','slightly higher than','Location','SouthEast');
grid on;
saveas(img5, sprintf('Tutorial_2_1_question_2c_1.png'));
%Lower Current: I_low doesn't produce spikes, the calculation and simulation confirm the minimum current threshold for spiking.
%Higher Current:I_high produces spikes, it verifies that currents above the threshold trigger spiking behavior.
Tau = C*R;
i_0 = [I_th*1.05 I_th*1.1 I_th*1.15 I_th*1.2 I_th*1.25 I_th*1.30 I_th*1.25 I_th*1.35 I_th*1.4 I_th*1.45];
frates = zeros(size(i_0));
k=1;
for I_0 = i_0
    
    [I_app_f,v_f]  = Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential);
   
    spike_count = 0;
    
    for n = 1:length(v_f)
        if v_f(n) == reset_potential;
            spike_count = spike_count + 1;
        end
    end
    rate = spike_count/2;
    frates(k) = rate;
    k=k+1;
end
img1=figure
scatter(i_0, frates, 'LineWidth', 2);
xlabel('I_app');
ylabel('Firing rate(Hz)');
title('Firing rate based on injected current 1c');
legend('1c','Location','SouthEast');
grid on;
saveas(img1, sprintf('Tutorial_2_1_question_1c.png'));

f=zeros(size(i_0));
j=1;
for I_0 = i_0
    if ((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - reset_potential) > 0) & ...
       ((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - threshold) > 0)
        
        firing_rate = 1 / (Tau * log((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - reset_potential) / ...
                       (leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - threshold)));
    
    end
    
    f(j) = firing_rate;
    j = j + 1;
end

img6=figure
scatter(i_0, frates, 'LineWidth', 2);
hold on
scatter(i_0, f, 'LineWidth', 2);
xlabel('I_app');
ylabel('Firing rate(Hz)');
title('Firing rate based on injected current');
legend('1c','1d','Location','SouthEast');
grid on;
saveas(img6, sprintf('Tutorial_2_1_question_2c_2.png'));
sigma=[0.1,0.2]

I_lower = I_th*0.99;
[I_app_lower_n,v_lower_n] = n_Iapp(I_lower, t,leak_potential,C,R,dt,threshold,reset_potential,sigma(1));
I_higher = I_th*1.01;
[I_app_higher_n,v_higher_n] = n_Iapp(I_higher, t,leak_potential,C,R,dt,threshold,reset_potential,sigma(1));

img7=figure;
plot(t(1:2001), v_lower(1:2001), 'b', 'LineWidth', 2);
hold on;
plot(t(1:2001), v_higher(1:2001), 'r', 'LineWidth', 2);
hold on
plot(t(1:2001), v_lower_n(1:2001), 'g', 'LineWidth', 2);
hold on;
plot(t(1:2001), v_higher_n(1:2001), 'c', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Membrane Potential (V)');
title('Membrane Potential with Different Applied Currents in presense of noise');
legend('slightly lower than','slightly higher than','slightly lower than-with noise','slightly higher than-with noise','Location','SouthEast');
grid on;
saveas(img7, sprintf('Tutorial_2_1_question_2c_3.png'));


n_frates = zeros(length(i_0),length(sigma));
k=1;
l=1;
for sigma_I=sigma
    
for I_0 = i_0
    
    [I_app_f,v_f]  = n_Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential,sigma_I);
   
    spike_count = 0;
    
    for n = 1:length(v_f)
        if v_f(n) == reset_potential;
            spike_count = spike_count + 1;
        end
    end
    rate = spike_count/2;
    n_frates(k,l) = rate;
    k=k+1;
end
k=1;
l=l+1;
end

f_n=zeros(length(i_0),length(sigma));
j=1;
w=1;
for sigma_I=sigma
for I_0 = i_0
    if ((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - reset_potential) > 0) & ...
       ((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - threshold) > 0)
        
        firing_rate = 1 / (Tau * log((leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - reset_potential) / ...
                       (leak_potential + R * Iapp(I_0, t, leak_potential, C, R, dt, threshold, reset_potential) - threshold)));
    
    end
    
    f_n(j,w) = firing_rate;
    j = j + 1;
end
j=1;
w=w+1
end

img8=figure
scatter(i_0, frates, 'SizeData', 52);
hold on
scatter(i_0, f, 'SizeData', 52);
hold on
scatter(i_0, n_frates(:,1), 'Marker', '*', 'SizeData', 52);
hold on
scatter(i_0, f_n(:,1), 'Marker', '*', 'SizeData', 52);
hold on
scatter(i_0, n_frates(:,2), 'Marker', 's', 'SizeData', 52);
hold on
scatter(i_0, f_n(:,2), 'Marker', 's', 'SizeData',52);
xlabel('I_app');
ylabel('Firing rate(Hz)');
title('Firing rate based on injected current');
legend('1c-sigma=0','1d-sigma=0','1c-sigma=0.1','1d-sigma=0.1','1c-sigma=0.2','1d-sigma=0.2','Location','SouthEast');
grid on;
saveas(img8, sprintf('Tutorial_2_1_question_2c_4.png'));