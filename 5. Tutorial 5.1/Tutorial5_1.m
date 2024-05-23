close all
clear
clc
%a
dt=0.1e-3;
time_vector=0:dt:4-dt;
%b,c
rate=[20,100,10,50];
% Generate presynaptic spike train as a Poisson process
nSpikes = zeros(size(time_vector));
 
for i = 1:length(rate)
  % Generate Poisson distributed spike counts for each interval
  nSpikes((i-1)*round(1/dt)+1:i*round(1/dt)) = (rand(1,round(1/dt))<rate(i)*dt);
end
f1=figure;
plot(time_vector,nSpikes)
xlabel('Time (Sec.)')
ylabel('Spikes')
title('Spikes-Part c');
grid on
saveas(f1, sprintf('Part_c.png'));
%d
 
dG=1e-9;
G_syn(1)=0;
for i=2:length(time_vector)
    G_syn(i)=G_syn(i-1)-dt*G_syn(i-1)/100e-3;
   if nSpikes(i)== 1;
       G_syn(i)=G_syn(i)+dG;
   end
end
f2=figure;
plot(time_vector,G_syn)
xlabel('Time (Sec.)')
ylabel('Conductance Vector')
title('Conductance Vector-Part d');
grid on
saveas(f2, sprintf('Part_d.png'));
%part e&f
p_0=0.5;
D(1)=1;
dG=1e-9;
G_syn_2(1)=0;
for i=2:length(time_vector)
   D(i)=D(i-1)+dt*(1-D(i-1))/0.25;
   G_syn_2(i)=G_syn_2(i-1)-dt*G_syn_2(i-1)/100e-3;
    if nSpikes(i)==1
       D(i)=D(i)-p_0*D(i);
       G_syn_2(i)=G_syn_2(i)+5e-9*p_0*D(i-1);
    end
end
f3=figure;
plot(time_vector,G_syn_2)
xlabel('Time (Sec.)')
ylabel('Conductance Vector')
title('Conductance Vector-Part f');
grid on
saveas(f3, sprintf('Part_f.png'));
%part G
p_0=0.2;
F(1)=1;
D(1)=1;
G_syn_3(1)=0;
for i=2:length(time_vector)
   F(i)=F(i-1)+dt*(1-F(i-1))/0.25;
   D(i)=D(i-1)+dt*(1-D(i-1))/0.25;
   G_syn_3(i)=G_syn_3(i-1)-dt*G_syn_3(i-1)/100e-3;
    if nSpikes(i)==1
      
       F(i)=F(i)+0.25*((1/p_0)-F(i-1));
       D(i)=D(i)-0.2*F(i-1)*D(i-1);
       G_syn_3(i)=G_syn_3(i)+4e-9*0.2*F(i-1)*D(i-1);
    end
end
f4=figure;
plot(time_vector,G_syn_3)
xlabel('Time (Sec.)')
ylabel('Conductance Vector')
title('Conductance Vector-Part h');
grid on
saveas(f4, sprintf('Part_h.png'));
