clc; clear all;
close all;

load('COVID_data.mat')  %load data file

peak1= 1:101;   %Feb 15 - May 25
peak2= 102:276; %May 26 - Nov 16

%modified Gaussian function
fun= @(p,x) p(1)+p(2).*exp(-((x-p(3)).^2)./(2.*(p(4).^2.*x)));

p01=[1 2100 30 100];    %for 1st peak of "currently infected" graph
pInfectedPeak1=lsqcurvefit(fun,p01,peak1,Australia.currently_infected(peak1));
p02=[110 3000 300 100]; %for 2nd peak of "currently infected" graph
pInfectedPeak2=lsqcurvefit(fun,p02,peak2,Australia.currently_infected(peak2));

p03=[0.1 2100 27 105];  %for 1st peak of "daily new cases" graph
pNewCasePeak1=lsqcurvefit(fun,p03,peak1,Australia.daily_new_cases(peak1));
p04= [110 4000 100 90]; %for 2nd peak of "daily new cases" graph
pNewCasePeak2=lsqcurvefit(fun,p04,peak2,Australia.daily_new_cases(peak2));


p05= [0.001 100 27 105];    %for 1st peak of "daily new death" graph
pNewDeathPeak1=lsqcurvefit(fun,p05,peak1,Australia.daily_new_death(peak1));
p06= [10 20 200 30];   %for 2nd peak of "daily new death" graph
pNewDeathPeak2=lsqcurvefit(fun,p06,peak2,Australia.daily_new_death(peak2));

subplot(1,3,1)  %place the 1st graph of 3
hold on
ylabel('Population');   % y axis
xlabel('Time (days)');  % x axis
title('Australia: Daily New Cases')
plot (Australia.daily_new_cases,'bo')   %plot daily new cases data with blue cirlces
plot (peak1, fun(pNewCasePeak1,peak1),'r-','LineWidth',2)   %plot with red lines for each peak
plot (peak2, fun(pNewCasePeak2,peak2),'r-','LineWidth',2)

subplot(1,3,2)  %place the 2nd graph of 3
hold on
ylabel('Population');   % y axis
ylim([0 70])
xlabel('Time (days)');  % x axis
title('Australia: Daily New Death')
plot (Australia.daily_new_death,'ko')   %plot daily new death data with black circles
plot (peak1, fun(pNewDeathPeak1,peak1),'r-','LineWidth',2)  %plot with red lines for each peak
plot (peak2, fun(pNewDeathPeak2,peak2),'r-','LineWidth',2)

subplot(1,3,3)  %place the 3rd graph of 3
hold on
ylabel('Population');   % y axis
xlabel('Time (days)');  % x axis
title('Australia: Currently Infected')  
plot (Australia.currently_infected,'go')    %plot currently infected data with green circles
plot (peak1, fun(pInfectedPeak1,peak1),'r-','LineWidth',2)  %plot with red lines for each peak
plot (peak2, fun(pInfectedPeak2,peak2),'r-','LineWidth',2)

%% using SIR model
clc; clear all;
close all;

load('COVID_data.mat')

N= 25499884; %australia population
I0= Australia.currently_infected(101); %initially infected on May 25th 
R0 = 1450000;   %initial recovered population        //change!
gamma = 0.3;   %duration of infectious period   //change!
r0 = 1.089;  %reproductive value              //change!
tspan = [101,276];  %length of time

recoveryRate = 1./gamma;   %recovery rate
beta = r0.*recoveryRate./N;   %infection rate
S0 = N -I0 -R0;   %initial susceptible population
N0 = I0 +R0 +S0;  %total population

pars = [beta, recoveryRate, N, r0]; %parameters
y0 = [S0 I0 R0];    %initial values

[t,y] = ode45(@SIR, tspan, y0, [], pars);   %SIR model

subplot(1,3,1)
hold on
ylabel('Population');   % y axis
xlabel('Time (days)');  % x axis
title('Australia: Susceptible')
plot(t,y(:,1),'b-','LineWidth',2)

subplot(1,3,2)
hold on
ylabel('Population');   % y axis
xlabel('Time (days)');  % x axis
title('Australia: Infectious')
plot(101:276,Australia.currently_infected(101:276),'go') %plot the given data(blue circle)
plot(t+30,y(:,2), 'r-','LineWidth',2);   %plot the infection fitting graph(red line)
legend('Currently Infected','SIR Infectious');

subplot(1,3,3)
hold on
ylabel('Population');   % y axis
xlabel('Time (days)');  % x axis
title('Australia: Recovered')
plot(t, y(:,3),'m-','LineWidth',2)

function f = SIR(t,y,pars)
f = zeros(3,1); %create three blanks
f(1) = -pars(1)*y(1)*y(2);  %ds/dt
f(2) = pars(1)*y(1)*y(2) -pars(2)*y(2);    %dI/dt
f(3) = pars(2)*y(2);  %dR/dt
end