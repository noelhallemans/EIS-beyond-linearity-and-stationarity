close all
clc
clear

%% Select experiment
temperature = 5;        % 5 or 25


if temperature == 5
    load('Data Samsung 48X\operandoEIS_5degrees.mat');    % Data measured at  5 degrees C
elseif temperature == 25
    load('Data Samsung 48X\operandoEIS_25degrees.mat');   % Data measured at 25 degrees C
end
disp("Temperature = " + temperature + " degrees C");

%% Variables in data file
% fs:   sampling frequency [Hz]
% Hexc: excited harmonics
% i:    applied odd random phase multisine current signals [A]
% P:    number of measured periods
% v:    measured voltage response [V]
% Zwt:  estimated time-varying impedance using operando EIS
% tZ:   time at which the impedance is computed

%% Derived parameters
N  = size(i,1);      % Number of measured data points per experiment
Ts = 1/fs;           % Sampling period [s]
T  = N*Ts;           % Measurment length [s]

Kexc = P*Hexc + 1;   % Indices of excited harmonics in spectra

t = 0:Ts:T-Ts;       % Time axis [s]
f = 0:fs/N:fs-fs/N;  % Frequency axis [Hz]

fmin = f(Kexc(1));   % Smalles excited frequency [Hz]



%% Show data in time domain

figure
subplot(211)
plot(t/60,i,'.');
ylabel('$i(t)$ [A]','Interpreter','latex','FontSize',11);  
xlim([0,T]/60); ylim([-2,7]);
title("Temperature $\approx$ " + temperature + " degrees C",'Interpreter','latex','FontSize',11);
subplot(212)
plot(t/60,v,'.','color',colors(2));
xlabel('time [min]','Interpreter','latex','FontSize',11);
ylabel('$v(t)$ [V]','Interpreter','latex','FontSize',11);  
xlim([0,T]/60); ylim([3,4.2])


%% Transform data into frequency domain
I = fft(i)/N;   % Discrete Fourier transform of the measured current
V = fft(v)/N;   % Discrete Fourier transform of the measured voltage

% Only retain frequencies under the Nyquist frequency
f = f(1:N/2); I = I(1:N/2); V = V(1:N/2);

Z = V(Kexc)./I(Kexc);   % Estimated impedance [Ohm]

fexc = f(Kexc);         % Excited frequencies

figure
stem(fexc,abs(I(Kexc)),'x','color',colors(1),'MarkerSize',8);
set(gca,'yscal','log'); set(gca,'xscal','log');
hold on
loglog(f,abs(I),'.','color',colors(1));
stem(fexc,abs(V(Kexc)),'.','color',colors(2),'MarkerSize',12);
set(gca,'yscal','log'); set(gca,'xscal','log');
hold on
loglog(f,abs(V),'.','color',colors(2));
xlim([1/T,100]); xticks(10.^(-4:2));
xlabel('frequency [Hz]','Interpreter','latex','FontSize',11);
ylabel('Magnitude [A,V]','Interpreter','latex','FontSize',11);
title("Temperature $\approx$ " + temperature + " degrees C",'Interpreter','latex','FontSize',11);

%% Time-varying impedance estimated with operando EIS

colormap = jet(length(tZ));    % Color map for time-varying impedance

figure
subplot(211)
box on
for ii = 50:10:length(tZ)-50
    loglog(fexc,1000*abs(Zwt(:,ii)),'. -','color',colormap(ii,:));
    hold on
end
if temperature == 5
    text(10,50,'$t=0$','Interpreter','latex','FontSize',9,'color',colormap(1,:));
    text(0.01,25,'$t=T$','Interpreter','latex','FontSize',9,'color',colormap(end,:));
elseif temperature == 25
    text(10,28,'$t=0$','Interpreter','latex','FontSize',9,'color',colormap(1,:));
    text(0.01,25,'$t=T$','Interpreter','latex','FontSize',9,'color',colormap(end,:));
end
xlim([fexc(1),100]);
ylabel('$\vert Z(\omega,t)\vert$ [m$\Omega$]','Interpreter','latex','FontSize',11);
title("Temperature $\approx$ " + temperature + " degrees C",'Interpreter','latex','FontSize',11);

subplot(212)
box on
for ii = 50:10:length(tZ)-50
    semilogx(fexc,180*angle(Zwt(:,ii))/pi,'. -','color',colormap(ii,:));
    hold on
end
xlim([fexc(1),100]);
xlabel('frequency [Hz]','Interpreter','latex','FontSize',11);
ylabel('$\angle Z(\omega,t)$ [$^\circ$]','Interpreter','latex','FontSize',11);

%% Nyquist plot
figure
box on
hold on
for ii = length(tZ)-50:-10:50
    plot(1000*real(Zwt(:,ii)),-1000*imag(Zwt(:,ii)),'. -','color',colormap(ii,:));  
end
axis equal
if temperature == 5
    text(10,50,'$t=0$','Interpreter','latex','FontSize',9,'color',colormap(1,:));
    text(0.01,25,'$t=T$','Interpreter','latex','FontSize',9,'color',colormap(end,:));
elseif temperature == 25
    text(10,28,'$t=0$','Interpreter','latex','FontSize',9,'color',colormap(1,:));
    text(0.01,25,'$t=T$','Interpreter','latex','FontSize',9,'color',colormap(end,:));
end
xlabel('$Z_\mathrm{r}(\omega,t)$ [m$\Omega$]','Interpreter','latex','FontSize',11);
ylabel('$-Z_\mathrm{j}(\omega,t)$ [m$\Omega$]','Interpreter','latex','FontSize',11);
title("Temperature $\approx$ " + temperature + " degrees C",'Interpreter','latex','FontSize',11);

%% Color function
function color = colors(i)
    MatlabColors = [0, 0.4470, 0.7410;
        0.8500, 0.3250, 0.0980];
    color = MatlabColors(i,:);
end

