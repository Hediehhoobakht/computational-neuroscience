clear
clc
% computational neuro science tutorial 3_2 mfile
VMAX = 50e-3; % V from tutorial 2_3 parameters
EL = -.07; % V
VTH = -.05; % V
VRESET = -.08; % V
DELTH = .002; % V
GL = 10e-9; % S
CM = 100e-12; % farad
A = 2e-9; % (Ohm)^-1
B = 0e-9; % Amp
TAUSRA = .15; % SRA period (s)
EK = -.08; % V
dt = .01e-3; % sampling time
T = 5e-3; % total time
% AELIF model 1st part
SIGMA = 50e-12;
SD = SIGMA/sqrt(dt);
t = 0:dt:100-dt; % initializing time vector
IAPP = randn(1,length(t))*SD; % initializing IAPP
V = zeros(size(IAPP)); % initiating j*n matrix for saving voltages
V(:,1) = EL; % setting initial condition of voltages
SPIKES = zeros(size(IAPP)); % initializing matrix of spikes
ISRA = zeros(size(IAPP)); % initializing I_SRA
for j = 1:size(IAPP,1)
    for i = 2:size(IAPP,2) % loop for forward integration through euler method
        V(j,i) = V(j,i-1)+dt*GL*(EL-V(j,i-1)+DELTH*exp((V(j,i-1)VTH)/DELTH))/CM+dt*(IAPP(j,i-1)-ISRA(j,i-1))/CM;


       ISRA(j,i)
= ISRA(j,i-1)+dt*(A*(V(j,i-1)-EL)-ISRA(j,i-1))/TAUSRA;


       if
(V(j,i)>VMAX)
%
if
potential
is
above
threshold


           SPIKES(j,i)
= 1; %
recording
spikes


           V(j,i)
= VRESET; %
potential
reset


           ISRA(j,i)
= ISRA(j,i) + B;




       end


   end

end

ISI
= diff(find(SPIKES))*dt*1e3; %
calculating
ISI

figure

histogram(ISI,25);
title('ISI
histogrm,
b
=
0
nA,
\mu
=
0
nA,
\sigma
=
50
pA');

xlabel('ISI
time
length
(ms)');
ylabel('N');

MISI
= mean(ISI); %
calculating
mean
ISI

SDISI
= std(ISI); %
calculating
standard
deviation
of
ISI

CV
= SDISI/MISI;
%
 calculating
CV

%

nhm
= 100e-3/dt;
%
time
step
width
for
100
ms



SPI_100 = sum(reshape(SPIKES,nhm,size(SPIKES,2)/nhm)); % number of spikes in each 100
ms
MSPI_100 = mean(SPI_100);
SDSPI_100 = std(SPI_100);
FNSPI_100 = SDSPI_100^2/MSPI_100; % fano factor
%
tw = 10e-3:10e-3:1; % initializing time windows
for i = 1:length(tw)
    nw = floor(tw(i)/dt); % time step width
    nnw = floor(length(t)/nw); % number of bins with the current time width
    SPI_tw = sum(reshape(SPIKES(1:nw*nnw),nw,nnw)); % number of spikes in each 100 ms
    MSPI(i) = mean(SPI_tw);
    SDSPI(i) = std(SPI_tw);
    FNSPI(i) = SDSPI(i)^2/MSPI(i); % fano factor
end
figure
plot(tw*1e3,MSPI); title('spike mean vs window width'); xlabel('time window width
(ms)'); ylabel('mean of spike');
figure
plot(tw*1e3,SDSPI); title('spike standard deviation vs window width'); xlabel('time
window width (ms)'); ylabel('SD of spike');
figure
plot(tw*1e3,FNSPI); title('Spike Fano factor vs window width, b = 0 nA, \mu = 0 nA,
\sigma = 50 pA'); xlabel('time window width (ms)'); ylabel('Fano factor of spike');
% 1-b
B = 1e-9; % Amp
IAPP = randn(1,length(t))*SD; % initializing IAPP
V = zeros(size(IAPP)); % initiating j*n matrix for saving voltages
V(:,1) = EL; % setting initial condition of voltages
SPIKES = zeros(size(IAPP)); % initializing matrix of spikes
ISRA = zeros(size(IAPP)); % initializing I_SRA
for j = 1:size(IAPP,1)
    for i = 2:size(IAPP,2) % loop for forward integration through euler method
        V(j,i) = V(j,i-1)+dt*GL*(EL-V(j,i-1)+DELTH*exp((V(j,i-1)VTH)/DELTH))/CM+dt*(IAPP(j,i-1)-ISRA(j,i-1))/CM;


       ISRA(j,i)
= ISRA(j,i-1)+dt*(A*(V(j,i-1)-EL)-ISRA(j,i-1))/TAUSRA;


       if
(V(j,i)>VMAX)
%
if
potential
is
above
threshold


           SPIKES(j,i)
= 1; %
recording
spikes


        
  V(j,i)
= VRESET; %
potential
reset


           ISRA(j,i)
= ISRA(j,i) + B;




       end


   end

end

ISI
= diff(find(SPIKES))*dt*1e3; %
calculating
ISI

figure

histogram(ISI,25);
title('ISI
histogrm,
b
=
1
nA,
\mu
=
0
nA,
\sigma
=
50
pA');

xlabel('ISI
time
(ms)');
ylabel('N');

MISI
= mean(ISI); %
calculating
mean
ISI



SDISI = std(ISI); % calculating standard deviation of ISI
CV = SDISI/MISI; %  calculating CV
%
nhm = 100e-3/dt; % time step width for 100 ms
SPI_100 = sum(reshape(SPIKES,nhm,size(SPIKES,2)/nhm)); % number of spikes in each 100
ms
MSPI_100 = mean(SPI_100);
SDSPI_100 = std(SPI_100);
FNSPI_100 = SDSPI_100^2/MSPI_100; % fano factor
%
tw = 10e-3:10e-3:1; % initializing time windows
for i = 1:length(tw)
    nw = floor(tw(i)/dt); % time step width
    nnw = floor(length(t)/nw); % number of bins with the current time width
    SPI_tw = sum(reshape(SPIKES(1:nw*nnw),nw,nnw)); % number of spikes in each 100 ms
    MSPI(i) = mean(SPI_tw);
    SDSPI(i) = std(SPI_tw);
    FNSPI(i) = SDSPI(i)^2/MSPI(i); % fano factor
end
figure
plot(tw*1e3,MSPI); title('spike mean vs window width'); xlabel('time window width
(ms)'); ylabel('mean of spike');
figure
plot(tw*1e3,SDSPI); title('spike standard deviation vs window width'); xlabel('time
window width (ms)'); ylabel('SD of spike');
figure
plot(tw*1e3,FNSPI); title('Spike Fano factor vs window width, b = 1 nA, \mu = 0 nA,
\sigma = 50 pA'); xlabel('time window width (ms)'); ylabel('Fano factor of spike');
% 1-c
B = 0e-9; % Amp
SIGMA = 20e-12;
SD = SIGMA/sqrt(dt);
t = 0:dt:100-dt; % initializing time vector
IAPP = randn(3,length(t))*SD+[0;.1;.2]*1e-9; % initializing IAPP
V = zeros(size(IAPP)); % initiating j*n matrix for saving voltages
V(:,1) = EL; % setting initial condition of voltages
SPIKES = zeros(size(IAPP)); % initializing matrix of spikes
ISRA = zeros(size(IAPP)); % initializing I_SRA
for j = 1:size(IAPP,1)
    for i = 2:size(IAPP,2) % loop for forward integration through euler method
        V(j,i) = V(j,i-1)+dt*GL*(EL-V(j,i-1)+DELTH*exp((V(j,i-1)VTH)/DELTH))/CM+dt*(IAPP(j,i-1)-ISRA(j,i-1))/CM;


       ISRA(j,i)
= ISRA(j,i-1)+dt*(A*(V(j,i-1)-EL)-ISRA(j,i-1))/TAUSRA;


       if
(V(j,i)>VMAX)
%
if
potential
is
above
threshold


           SPIKES(j,i)
= 1; %
recording
spikes


           V(j,i)
= VRESET; %
potential
reset


           ISRA(j,i)
= ISRA(j,i) + B;




       end



    end
end
ISI0 = diff(find(SPIKES(1,:)))*dt*1e3; % calculating ISI
ISI1 = diff(find(SPIKES(2,:)))*dt*1e3; % calculating ISI
ISI2 = diff(find(SPIKES(3,:)))*dt*1e3; % calculating ISI
figure
histogram(ISI0,25); title('ISI histogrm, b = 0 nA, \mu = 0 nA, \sigma = 20 pA');
xlabel('ISI time (ms)'); ylabel('N');
figure
histogram(ISI1,25); title('ISI histogrm, b = 0 nA, \mu = 0.1 nA, \sigma = 20 pA');
xlabel('ISI time (ms)'); ylabel('N');
figure
histogram(ISI2,25); title('ISI histogrm, b = 0 nA, \mu = 0.2 nA, \sigma = 20 pA');
xlabel('ISI time (ms)'); ylabel('N');
MISI0 = mean(ISI0); % calculating mean ISI
SDISI0 = std(ISI0); % calculating standard deviation of ISI
CV0 = SDISI0/MISI0; %  calculating CV
MISI1 = mean(ISI1); % calculating mean ISI
SDISI1 = std(ISI1); % calculating standard deviation of ISI
CV1 = SDISI1/MISI1; %  calculating CV
MISI2 = mean(ISI2); % calculating mean ISI
SDISI2 = std(ISI2); % calculating standard deviation of ISI
CV2 = SDISI2/MISI2; %  calculating CV
%
nhm = 100e-3/dt; % time step width for 100 ms
SPI0_100 = sum(reshape(SPIKES(1,:),nhm,size(SPIKES(1,:),2)/nhm)); % number of spikes
in each 100 ms
MSPI0_100 = mean(SPI0_100);
SDSPI0_100 = std(SPI0_100);
FNSPI0_100 = SDSPI0_100^2/MSPI0_100; % fano factor for mean = 0
nhm = 100e-3/dt; % time step width for 100 ms
SPI1_100 = sum(reshape(SPIKES(2,:),nhm,size(SPIKES(2,:),2)/nhm)); % number of spikes
in each 100 ms
MSPI1_100 = mean(SPI1_100);
SDSPI1_100 = std(SPI1_100);
FNSP1_100 = SDSPI1_100^2/MSPI1_100; % fano factor for mean = .1 nA
nhm = 100e-3/dt; % time step width for 100 ms
SPI2_100 = sum(reshape(SPIKES(3,:),nhm,size(SPIKES(3,:),2)/nhm)); % number of spikes
in each 100 ms
MSPI2_100 = mean(SPI2_100);
SDSPI2_100 = std(SPI2_100);
FNSPI2_100 = SDSPI2_100^2/MSPI2_100; % fano factor for mean = .2 nA
% part B
dt = .1e-3; % s
rate = 20; % Hz
t = 0:dt:100-dt; % initializing time
SPIKES = rand(1,length(t))<rate*dt;
ISI = diff(find(SPIKES))*dt*1e3; % calculating ISI
figure 
histogram(ISI,25); title('ISI histogrm from poisson distribution, rate = 20 Hz,
\Delta{t} = 0.1 ms'); xlabel('ISI time length (ms)'); ylabel('N');
MISI = mean(ISI); % calculating mean ISI
SDISI = std(ISI); % calculating standard deviation of ISI
CV = SDISI/MISI; %  calculating CV
%
nhm = 100e-3/dt; % time step width for 100 ms
SPI_100 = sum(reshape(SPIKES,nhm,size(SPIKES,2)/nhm)); % number of spikes in each 100
ms
MSPI_100 = mean(SPI_100);
SDSPI_100 = std(SPI_100);
FNSPI_100 = SDSPI_100^2/MSPI_100; % fano factor
%
tw = 10e-3:10e-3:1; % initializing time windows
for i = 1:length(tw)
    nw = floor(tw(i)/dt); % time step width
    nnw = floor(length(t)/nw); % number of bins with the current time width
    SPI_tw = sum(reshape(SPIKES(1:nw*nnw),nw,nnw)); % number of spikes in each 100 ms
    MSPI(i) = mean(SPI_tw);
    SDSPI(i) = std(SPI_tw);
    FNSPI(i) = SDSPI(i)^2/MSPI(i); % fano factor
end
figure
plot(tw*1e3,FNSPI); title({'Spike Fano factor vs window width','from poisson
distribution, rate = 20 Hz, \Delta{t} = 0.1 ms'}); xlabel('time window width (ms)');
ylabel('Fano factor of spike');
%
t = 0:dt:10-dt; % initializing time
SPIKES = rand(1000,length(t))<rate*dt; % producing spikes
CSPIKES = cumsum(SPIKES,1);
MCSPIKES = mean(CSPIKES,1); % calculating row mean
VCSPIKES = var(CSPIKES,0,1); % calculating row variance
FNCSPIKES = VCSPIKES./MCSPIKES; % calculatinf Fano factor
figure
plot(t*1e3,FNCSPIKES); title({'Spike Fano factor of 1000 trials vs elapsed time','from
poisson distribution, rate = 20 Hz, \Delta{t} = 0.1 ms'}); xlabel('elapsed time
(ms)'); ylabel('Fano factor of spike');
%
nhm = 200e-3/dt; % time step width for 100 ms
for i = 1:1000
SPI_200(i,:) = sum(reshape(SPIKES(i,:),nhm,size(SPIKES(i,:),2)/nhm)); % number of
spikes in each 100 ms
end
MSPI_200 = mean(SPI_200,1); % mean of each column
VSPI_200 = var(SPI_200,0,1); % var of each column
FNSPI_200 = VSPI_200./MSPI_200; % calculatinf Fano factor
figure
plot(FNSPI_200); title({'Spike Fano factor of 1000 trials vs elapsed "200 ms time
bins"','from poisson distribution, rate = 20 Hz, \Delta{t} = 0.1 ms'}); xlabel('number 
of elapsed 200 ms time bins'); ylabel('Fano factor of spike'); 
%
for i = 1:length(t)
    VSPIKES(:,i) = rand(1000,1)<(25+20*sin(2*pi*t(i)))*dt; % producing spikes
end
CVSPIKES = cumsum(VSPIKES,1);
MCVSPIKES = mean(CVSPIKES,1); % calculating row mean
VCVSPIKES = var(CVSPIKES,0,1); % calculating row variance
FNCVSPIKES = VCVSPIKES./MCVSPIKES; % calculatinf Fano factor
figure
plot(t*1e3,FNCVSPIKES); title({'Spike Fano factor of 1000 trials vs elapsed
time','from poisson distribution with variable rate'}); xlabel('elapsed time (ms)');
ylabel('Fano factor of spike');
%
nhm = 200e-3/dt; % time step width for 100 ms
for i = 1:1000
VSPI_200(i,:) = sum(reshape(VSPIKES(i,:),nhm,size(VSPIKES(i,:),2)/nhm)); % number of
spikes in each 100 ms
end
MVSPI_200 = mean(VSPI_200,1); % mean of each column
VVSPI_200 = var(VSPI_200,0,1); % var of each column
FNVSPI_200 = VVSPI_200./MSPI_200; % calculatinf Fano factor
figure
plot(FNVSPI_200); title({'Spike Fano factor of 1000 trials vs elapsed "200 ms time bins','from poisson distribution with variable rate'}); 
xlabel('number of elapsed 200ms time bins'); ylabel('Fano factor of spike'); 