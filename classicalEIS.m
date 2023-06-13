close all
clc
clear

%% Classical EIS experiments on Samsung 48X cells at different SOC

%% Select experiment
temperature = 5;        % 5 or 25
SOC_index   = 1;        
% For temperature =  5: 1=10%, 2=20%, ..., 8=80% 
% For temperature = 25: 1=10%, 2=20%, ..., 9=90% 

if temperature == 5
    load('Data Samsung 48X\classicalEIS_5degrees.mat');    % Data measured at  5 degrees C
elseif temperature == 25
    load('Data Samsung 48X\classicalEIS_25degrees.mat');   % Data measured at 25 degrees C
end
SOC = SOC(SOC_index);          
disp("Temperature = " + temperature + " degrees C,  SOC = " + SOC + "%");

%% Variables in data file
% fs:   sampling frequency [Hz]
% Hexc: excited harmonics
% i:    applied odd random phase multisine current signals [A]
% P:    number of measured periods
% SOC:  SOC levels at which data has been measured [%]
% v:    measured voltage response [V]

i = i(:,SOC_index);  % Current signal at selected SOC [A]
v = v(:,SOC_index);  % Voltage signal at selected SOC [V]


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
xlim([0,T]/60); ylim([-3,3]);
title("Temperature = " + temperature + " degrees C,  SOC = " + SOC + "\%",'Interpreter','latex','FontSize',11);
subplot(212)
plot(t/60,v,'.','color',[0.8500 0.3250 0.0980]);
xlabel('time [min]','Interpreter','latex','FontSize',11);
ylabel('$v(t)$ [V]','Interpreter','latex','FontSize',11);  
xlim([0,T]/60); ylim([3,4.2])


%% Transform data into frequency domain
I = fft(i)/N;   % Discrete Fourier transform of the measured current
V = fft(v)/N;   % Discrete Fourier transform of the measured voltage

% Only retain frequencies under the Nyquist frequency
f = f(1:N/2); I = I(1:N/2); V = V(1:N/2);

Z = V(Kexc)./I(Kexc);   % Estimated impedance [Ohm]

fexc = f(Kexc);         % Excited frequencies [Hz]

figure
subplot(4,1,[1,2])
stem(fexc,abs(I(Kexc)),'x','color',colors(1),'MarkerSize',8);
set(gca,'yscal','log'); set(gca,'xscal','log');
hold on
loglog(f,abs(I),'.','color',colors(1));
stem(fexc,abs(V(Kexc)),'.','color',colors(2),'MarkerSize',12);
set(gca,'yscal','log'); set(gca,'xscal','log');
hold on
loglog(f,abs(V),'.','color',colors(2));
text(0.001,0.2,'$I(k)$','color',colors(1),'Interpreter','latex','FontSize',11);
text(0.001,0.008,'$V(k)$','color',colors(2),'Interpreter','latex','FontSize',11);
xlim([1/T,100]); 
ylabel('Magnitude [A,V]','Interpreter','latex','FontSize',11);
title("Temperature = " + temperature + " degrees C,  SOC = " + SOC + "\%",'Interpreter','latex','FontSize',11);

subplot(413);
loglog(fexc,1000*abs(Z),'. -k','MarkerSize',10);
xlim([1/T,100]); 
ylabel('$\vert Z(\omega_k)\vert$ [m$\Omega$]','Interpreter','latex','FontSize',11);  
subplot(414);
semilogx(fexc,180*angle(Z)/pi,'. -k','MarkerSize',10);
xlim([1/T,100]); 
xlabel('frequency [Hz]','Interpreter','latex','FontSize',11);
ylabel('$\angle Z(\omega_k)\vert$ [$^\circ$]','Interpreter','latex','FontSize',11);  


figure
plot(1000*real(Z),-1000*imag(Z),'. -k','MarkerSize',10);
axis equal
xlabel('$Z_\mathrm{r}(\omega)$ [m$\Omega$]','Interpreter','latex','FontSize',11);
ylabel('$-Z_\mathrm{j}(\omega)$ [m$\Omega$]','Interpreter','latex','FontSize',11); 
title("Temperature = " + temperature + " degrees C,  SOC = " + SOC + "\%",'Interpreter','latex','FontSize',11);


%% Color function
function color = colors(i)
    MatlabColors = [0, 0.4470, 0.7410;
        0.8500, 0.3250, 0.0980];
    color = MatlabColors(i,:);
end
