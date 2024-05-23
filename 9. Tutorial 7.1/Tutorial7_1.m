close all
clear
clc

r_max=100;
theta_E=-5;
theta_I=0;
alpha_E=0.05;
alpha_I=1;
W_EE=2;
W_EI=2.5;
W_IE=-2.5;
W_II=-2;
t_max=3;
r_E(1)=50;
r_I(1)=50;
dt=0.1e-3;
t=0:dt:t_max;
Istim = zeros(size(t));
ind = find (t > 1 & t < 2);
Istim(ind) = 20;
Ibase_E = 0;
Ibase_I = 0;
tau_E=5e-3;
tau_I=5e-3;
Iapp_E = ones(size(t))*Ibase_E;
Iapp_I = ones(size(t))*Ibase_I+Istim;
r2null = ones(size(t))
for n = 2:length(t)
    I_E(n)=W_EE*r_E(n-1)+W_IE*r_I(n-1)+Iapp_E(n-1);
    I_I(n)=W_EI*r_E(n-1)+W_II*r_I(n-1)+Iapp_I(n-1);
    r_E_temp(n) =  r_E(n-1)+(dt/tau_E)*(-r_E(n-1)+alpha_E*((I_E(n)-theta_E)^2).*sign(I_E(n)-theta_E));
    r_I_temp(n) =  r_I(n-1)+(dt/tau_I)*(-r_I(n-1)+alpha_I*((I_I(n)-theta_I)));
    r_E(n) = min(max(r_E_temp(n), 0), r_max);
    r_I(n) = min(max(r_I_temp(n), 0), r_max);
    %r2null(n) = (1/tau_E)*(-r_E(n-1)+alpha_E*((I_E(n)-theta_E)^2).*sign(I_E(n)-theta_E));

end
% allowed_vals = find((r_E >= 0 ).*( r_E <= r_max ));
% allowed_vals1 = find((r_I>= 0 ).*( r_I <= r_max ));

f1=figure(1);
% subplot(3,1,1);
% plot(t,r_E);
% xlabel('t (sec)')
% ylabel('r_{E} (Hz)')
% grid on
% hold on
% plot(t,r_I); 
% 
% xlabel('t (sec)')
% ylabel('r_{I} (Hz)')
% plot(t(allowed_vals),r_E(allowed_vals));
subplot(1,2,1)
plot(t,r_E);
xlabel('t (sec)')
ylabel('r (Hz)')
grid on
hold on
%plot(t(allowed_vals1),r_I(allowed_vals1)); 
plot(t,r_I); 
xlabel('t (sec)')
ylabel('r (Hz)')
legend('r_{E}','r_{I}')
grid on

subplot(1,2,2)
plot(r_I,r_E)
ylabel('r_{E}')
xlabel('r_{I}')
% allowed_vals2 = find((r2null > 0 ).*( r2null < r_max ));
% figure(2)
% plot(r_E(allowed_vals2),r2null(allowed_vals2),'k');
% subplot(3,1,3);
% r_E1=1:length(r_I(allowed_vals1))
% plot(t(allowed_vals1),r_I(allowed_vals1))
% hold on 
% r_I1=1:length(r_I(allowed_vals))
% plot(r_E(allowed_vals),t(allowed_vals))
% axis([0 t_max 0 r_max])
% xlabel('r_{E} (sec)')
% ylabel('r_{I} (Hz)')
grid on
suptitle('part 1');

saveas(f1, sprintf('1.png'));

%% part2
r_max=100;
theta_E=-5;
theta_I=0;
alpha_E=0.05;
alpha_I=1;
W_EE=2;
W_EI=2.5;
W_IE=-2.5;
W_II=-2;
t_max=3;
r_E(1)=50;
r_I(1)=50;
dt=0.1e-3;
t=0:dt:t_max;
Istim = zeros(size(t));
ind = find (t > 1 & t < 2);
Istim(ind) = 20;
Ibase_E = 25;
Ibase_I = 15;
tau_E=5e-3;
tau_I=5e-3;
Iapp_E = ones(size(t))*Ibase_E;
Iapp_I = ones(size(t))*Ibase_I+Istim;

for n = 2:length(t)
    I_E(n)=W_EE*r_E(n-1)+W_IE*r_I(n-1)+Iapp_E(n-1);
    I_I(n)=W_EI*r_E(n-1)+W_II*r_I(n-1)+Iapp_I(n-1);
    r_E_temp(n) =  r_E(n-1)+(dt/tau_E)*(-r_E(n-1)+alpha_E*((I_E(n)-theta_E)^2).*sign(I_E(n)-theta_E));
    r_I_temp(n) =  r_I(n-1)+(dt/tau_I)*(-r_I(n-1)+alpha_I*((I_I(n)-theta_I)));
    r_E(n) = min(max(r_E_temp(n), 0), r_max);
    r_I(n) = min(max(r_I_temp(n), 0), r_max);
    
end
allowed_vals = find((r_E >= 0 ).*( r_E <= r_max ));
allowed_vals1 = find((r_I>= 0 ).*( r_I <= r_max ));

f2=figure(2);

subplot(1,2,1)
plot(t,r_E);
xlabel('t (sec)')
ylabel('r (Hz)')
grid on
hold on
%plot(t(allowed_vals1),r_I(allowed_vals1)); 
plot(t,r_I); 
xlabel('t (sec)')
ylabel('r (Hz)')
legend('r_{E}','r_{I}')
grid on

subplot(1,2,2)
plot(r_I,r_E)
ylabel('r_{E}')
xlabel('r_{I}')
% subplot(3,1,3);
% r_E1=1:length(r_I(allowed_vals1))
% plot(t(allowed_vals1),r_I(allowed_vals1))
% hold on 
% r_I1=1:length(r_I(allowed_vals))
% plot(r_E(allowed_vals),t(allowed_vals))
% axis([0 t_max 0 r_max])
% xlabel('r_{E} (sec)')
% ylabel('r_{I} (Hz)')
grid on
title('part 2');
%legend('r_{E}','r_{I}')
saveas(f2, sprintf('2.png'));
%% part 3
r_max=100;
theta_E=-5;
theta_I=0;
alpha_E=0.05;
alpha_I=1;
W_EE=2;
W_EI=2.5;
W_IE=-2.5;
W_II=-2;
t_max=3;
r_E(1)=50;
r_I(1)=50;
dt=0.1e-3;
t=0:dt:t_max;
Istim = zeros(size(t));
ind = find (t > 1 & t < 2);
Istim(ind) = 20;
Ibase_E = 0;
Ibase_I = 0;
tau_E=2e-3;
tau_I=10e-3;
Iapp_E = ones(size(t))*Ibase_E;
Iapp_I = ones(size(t))*Ibase_I+Istim;

for n = 2:length(t)
    I_E(n)=W_EE*r_E(n-1)+W_IE*r_I(n-1)+Iapp_E(n-1);
    I_I(n)=W_EI*r_E(n-1)+W_II*r_I(n-1)+Iapp_I(n-1);
    r_E_temp(n) =  r_E(n-1)+(dt/tau_E)*(-r_E(n-1)+alpha_E*((I_E(n)-theta_E)^2).*sign(I_E(n)-theta_E));
    r_I_temp(n) =  r_I(n-1)+(dt/tau_I)*(-r_I(n-1)+alpha_I*((I_I(n)-theta_I)));
    r_E(n) = min(max(r_E_temp(n), 0), r_max);
    r_I(n) = min(max(r_I_temp(n), 0), r_max);
    
end
allowed_vals = find((r_E >= 0 ).*( r_E <= r_max ));
allowed_vals1 = find((r_I>= 0 ).*( r_I <= r_max ));

f3=figure(3);

subplot(1,2,1)
plot(t,r_E);
xlabel('t (sec)')
ylabel('r (Hz)')
grid on
hold on
%plot(t(allowed_vals1),r_I(allowed_vals1)); 
plot(t,r_I); 
xlabel('t (sec)')
ylabel('r (Hz)')
legend('r_{E}','r_{I}')
grid on

subplot(1,2,2)
plot(r_I,r_E)
ylabel('r_{E}')
xlabel('r_{I}')
% subplot(3,1,3);
% r_E1=1:length(r_I(allowed_vals1))
% plot(t(allowed_vals1),r_I(allowed_vals1))
% hold on 
% r_I1=1:length(r_I(allowed_vals))
% plot(r_E(allowed_vals),t(allowed_vals))
% axis([0 t_max 0 r_max])
% xlabel('r_{E} (sec)')
% ylabel('r_{I} (Hz)')
grid on
title('part 3');
%legend('r_{E}','r_{I}')
saveas(f3, sprintf('3.png'));
%% part 4
r_max=100;
theta_E=-5;
theta_I=0;
alpha_E=0.05;
alpha_I=1;
W_EE=2;
W_EI=2.5;
W_IE=-2.5;
W_II=-2;
t_max=3;
r_E(1)=50;
r_I(1)=50;
dt=0.1e-3;
t=0:dt:t_max;
Istim = zeros(size(t));
ind = find (t > 1 & t < 2);
Istim(ind) = 20;
Ibase_E = 25;
Ibase_I = 15;
tau_E=2e-3;
tau_I=10e-3;
Iapp_E = ones(size(t))*Ibase_E;
Iapp_I = ones(size(t))*Ibase_I+Istim;

for n = 2:length(t)
    I_E(n)=W_EE*r_E(n-1)+W_IE*r_I(n-1)+Iapp_E(n-1);
    I_I(n)=W_EI*r_E(n-1)+W_II*r_I(n-1)+Iapp_I(n-1);
    r_E_temp(n) =  r_E(n-1)+(dt/tau_E)*(-r_E(n-1)+alpha_E*((I_E(n)-theta_E)^2).*sign(I_E(n)-theta_E));
    r_I_temp(n) =  r_I(n-1)+(dt/tau_I)*(-r_I(n-1)+alpha_I*((I_I(n)-theta_I)));
    r_E(n) = min(max(r_E_temp(n), 0), r_max);
    r_I(n) = min(max(r_I_temp(n), 0), r_max);
    
end
allowed_vals = find((r_E >= 0 ).*( r_E <= r_max ));
allowed_vals1 = find((r_I>= 0 ).*( r_I <= r_max ));

f4=figure(4);
subplot(1,2,1)
plot(t,r_E);
xlabel('t (sec)')
ylabel('r (Hz)')
grid on
hold on
%plot(t(allowed_vals1),r_I(allowed_vals1)); 
plot(t,r_I); 
xlabel('t (sec)')
ylabel('r (Hz)')
legend('r_{E}','r_{I}')
grid on

subplot(1,2,2)
plot(r_I,r_E)
ylabel('r_{E}')
xlabel('r_{I}')
% subplot(3,1,3);
% plot(t(allowed_vals),r_E(allowed_vals));
% hold on
% plot(t(allowed_vals1),r_I(allowed_vals1));

grid on
title('part 4');
%legend('r_{E}','r_{I}')
saveas(f4, sprintf('4.png'));