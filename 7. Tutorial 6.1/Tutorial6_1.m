close all
clear
clc
%part 1
D=1;
alpha_0 = 0.5;
WEE = 8;
p_r = 1;
tau_s = 2E-3;
x = 1.2;

r_0 = 0.1;
r_max = 100;
sig = 0.5;
tau_r = 10E-3;
%%part 1a
n = 1;
for s = 0:0.001:1/WEE
   s1(n) = s;
   S = WEE*s;
   f_S(n) = f(S,r_0,r_max,sig,x);
   n = n+1;
end
f1=figure(1)
hold off
plot(s1, f_S)

% xlabel('s')
% ylabel('f(WEE*S)')
% title('part 1.a')
% saveas(f1, sprintf('1A.png'));
%now plot s(r)
n = 1;
clear s;
for r = 0:0.01:r_max
   r1(n) = r;
   D = 1;
   s(n) = alpha_0*D*p_r*r*tau_s/(1+alpha_0*D*p_r*r*tau_s);
   n = n+1;
end
hold on
plot(s,r1,'r--');
xlabel('s')
ylabel('r (Hz)')
xlim([0 0.03])
legend('r = f(S)', 's(r), synaptic strength', 'location', 'northwest')
title('part 1.a')
saveas(f1, sprintf('1A.png'));
%part b
dt=0.1e-3;
t=0:dt:20;
clear s;clear r;
s(1)=0;r(1)=0;
for i=2:length(t)
    S=WEE*s(i-1);
    r(i)=r(i-1)+dt*(-r(i-1)+f(S,r_0,r_max,sig,x))/tau_r;
    s(i)=s(i-1)+dt*(-s(i-1)/tau_s+alpha_0*D*p_r*r(i)*(1-s(i-1)));
    if t(i)>=10 && t(i)<=10.05
        s(i)=0.05;
    end
end
f2=figure(2)
plot(t,r)
xlabel('t(sec)')
ylabel('r(t)(Hz)')
title('part 1-b')
saveas(f2, sprintf('1B.png'));
f3=figure(3)
plot(t,s)
xlabel('t(sec)')
ylabel('s(t)')
title('part 1-b')
saveas(f3, sprintf('1B1.png'));
%% part 2
tau_D=250e-3;
p_r=0.2;
r_0=0.1;
alpha_0=0.5;
WEE=60;
clear s; clear r;
n = 1;
for s = 0:0.001:1/WEE
   s1(n) = s;
   S = WEE*s;
   f_S(n) = f(S,r_0,r_max,sig,x);
   n = n+1;
end
f4=figure(4)
hold off
plot(s1, f_S)
n = 1;
clear s;clear r;
for r = 0:0.01:r_max
   r1(n) = r;
   D = 1/(1+p_r*r*tau_D);
   s(n) = alpha_0*D*p_r*r*tau_s/(1+alpha_0*D*p_r*r*tau_s);
   n = n+1;
end
hold on
plot(s,r1,'r--');
xlim([0 0.003])
xlabel('s')
ylabel('r (Hz)')
legend('r = f(S)', 's(r), synaptic strength', 'location', 'northwest')
title('part 2-a')
saveas(f4, sprintf('2A.png'));
%part 2-b
dt=0.1e-3;
t=0:dt:20;
clear s;clear r;
s(1)=0;r(1)=0;D(1)=1;
for i=2:length(t)
    S=WEE*s(i-1);
    r(i)=r(i-1)+dt*(-r(i-1)+f(S,r_0,r_max,sig,x))/tau_r;
    D(i) =D(i-1)+dt*(((1-D(i-1))/tau_D)-p_r*D(i-1)*r(i));
    s(i)=s(i-1)+dt*(-s(i-1)/tau_s+alpha_0*D(i)*p_r*r(i)*(1-s(i-1)));
    if t(i)>=10 && t(i)<=12
        s(i)=0.002;
    end
end
f5=figure(5)
plot(t,r)
xlabel('t(sec)')
ylabel('r(t)(Hz)')
title('part 2-b')
saveas(f5, sprintf('2B.png'));
f6=figure(6)
plot(t,s)
xlabel('t(sec)')
ylabel('s(t)')
title('part 2-b')
saveas(f6, sprintf('2B1.png'));
%% part 3
p_r=0.5;
WEE=35;
r_0=-0.1;
clear s; clear r;
n = 1;
for s = 0:0.001:1/WEE
   s1(n) = s;
   S = WEE*s;
   f_S(n) = f(S,r_0,r_max,sig,x);
   n = n+1;
end
f7=figure(7)
hold off
plot(s1, f_S)
n = 1;
clear s;clear r;
for r = 0:0.01:r_max
   r1(n) = r;
   D = 1/(1+p_r*r*tau_D);
   s(n) = alpha_0*D*p_r*r*tau_s/(1+alpha_0*D*p_r*r*tau_s);
   n = n+1;
end
hold on
plot(s,r1,'r--');
xlim([0 0.003])
xlabel('s')
ylabel('r (Hz)')
legend('r = f(S)', 's(r), synaptic strength', 'location', 'northwest')
title('part 3-a')
saveas(f7, sprintf('3A.png'));
dt=0.1e-3;
t=0:dt:20;
clear s;clear r;
s(1)=0;r(1)=0;D(1)=1;
for i=2:length(t)
    S=WEE*s(i-1);
    r(i)=r(i-1)+dt*(-r(i-1)+f(S,r_0,r_max,sig,x))/tau_r;
    D(i) =D(i-1)+dt*(((1-D(i-1))/tau_D)-p_r*D(i-1)*r(i));
    s(i)=s(i-1)+dt*(-s(i-1)/tau_s+alpha_0*D(i)*p_r*r(i)*(1-s(i-1)));
    if t(i)>=10 && t(i)<=10.6
        s(i)=0.002;
    end
end
f8=figure(8)
plot(t,r)
xlabel('t(sec)')
ylabel('r(t)(Hz)')
title('part 3-b')
saveas(f8, sprintf('3B.png'));
f9=figure(9)
plot(t,s)
xlabel('t(sec)')
ylabel('s(t)')
title('part 3-b')
saveas(f9, sprintf('3B1.png'));
%%part 3-c

dt=0.1e-3;
t=0:dt:20;
clear s;clear r;

r(1)=9;D(1)=1/(1+p_r*r(1)*tau_D);s(1)=alpha_0*D(1)*p_r*r(1)*tau_s/(1+alpha_0*D(1)*p_r*r(1)*tau_s);
for i=2:length(t)
    S=WEE*s(i-1);
    r(i)=r(i-1)+dt*(-r(i-1)+f(S,r_0,r_max,sig,x))/tau_r;
    D(i) =D(i-1)+dt*(((1-D(i-1))/tau_D)-p_r*D(i-1)*r(i));
    s(i)=s(i-1)+dt*(-s(i-1)/tau_s+alpha_0*D(i)*p_r*r(i)*(1-s(i-1)));
    if t(i)>=10 && t(i)<=10.6
        s(i)=0.002;
    end
end
f10=figure(10)
plot(t,r)
xlabel('t(sec)')
ylabel('r(t)(Hz)')
title('part 3-c')
saveas(f10, sprintf('3C.png'));
f11=figure(11)
plot(t,s)
xlabel('t(sec)')
ylabel('s(t)')
title('part 3-c')
saveas(f11, sprintf('3C1.png'));
%% part 4
tau_D=0.125
alpha_0=0.25;
p_r=1;
clear s; clear r;
r(1)=9;D(1)=1/(1+p_r*r(1)*tau_D);s(1)=alpha_0*D(1)*p_r*r(1)*tau_s/(1+alpha_0*D(1)*p_r*r(1)*tau_s);
n = 1;
for s = 0:0.001:1/WEE
   s1(n) = s;
   S = WEE*s;
   f_S(n) = f(S,r_0,r_max,sig,x);
   n = n+1;
end
f12=figure(12);
hold off
plot(s1, f_S)
n = 1;
clear s;clear r;
for r = 0:0.01:r_max
   r1(n) = r;
   D = 1/(1+p_r*r*tau_D);
   s(n) = alpha_0*D*p_r*r*tau_s/(1+alpha_0*D*p_r*r*tau_s);
   n = n+1;
end
hold on
plot(s,r1,'r--');
xlim([0 0.003])
xlabel('s')
ylabel('r (Hz)')
legend('r = f(S)', 's(r), synaptic strength', 'location', 'northwest')
title('part 4-a')
saveas(f12, sprintf('4A.png'));
dt=0.1e-3;
t=0:dt:20;
clear s;clear r;
s(1)=0;r(1)=0;D(1)=1;
for i=2:length(t)
    S=WEE*s(i-1);
    r(i)=r(i-1)+dt*(-r(i-1)+f(S,r_0,r_max,sig,x))/tau_r;
    D(i) =D(i-1)+dt*(((1-D(i-1))/tau_D)-p_r*D(i-1)*r(i));
    s(i)=s(i-1)+dt*(-s(i-1)/tau_s+alpha_0*D(i)*p_r*r(i)*(1-s(i-1)));
    if t(i)>=10 && t(i)<=10.6
        s(i)=0.002;
    end
end
f13=figure(13);
plot(t,r)
xlabel('t(sec)')
ylabel('r(t)(Hz)')
title('part 4-b')
saveas(f13, sprintf('4B.png'));
f14=figure(14)
plot(t,s)
xlabel('t(sec)')
ylabel('s(t)')
title('part 4-b')
saveas(f14, sprintf('4B1.png'));