close all
clc
clear


%% User defined parameters
fs    = 200;               % Sampling frequency [Hz]
fmax  = 80;                % Largest excited frequency [Hz], should be smaller than fs/2  
Tp    = 3*60;              % Period length [s]
irms  = 0.8;               % RMS value of the multisine [A]

decay = 0;                 % Decay of the multisine amplitudes over frequency (0 =  no decay)
L             = 4;         % One out of L consecutive odd harmonics will be removed
NexcperDecade = 22;        % Number of excited frequencies per decade


%% Derived and fixed parameters
fmin = 1/Tp;               % Smallest excited frequency [Hz]
Ts   = 1/fs;               % Sampling period [s]
Npp  = Tp/Ts;              % Number of points per period

t = 0:Ts:(Npp-1)*Ts;       % Time axis [s]
f = 0:fs/Npp:fs-fs/Npp;    % Frequency axis [Hz]

%% Generate vector with logarithmically distributed odd harmonics Hexc 
hexcmax = floor(fmax/fmin);   % Highest excited harmonic
a = hexcmax^(1/hexcmax);
Ndecades = log10(fmax/fmin);  

Hexc = linspace(0,hexcmax,NexcperDecade*Ndecades);
Hexc = round((a.^Hexc)/2)*2-1; 
Hexc = unique(Hexc);
Hexc(Hexc < 1) = [];

% Remove 1 out of L consecutive lines 
Ntmp = mod(length(Hexc),L);
Hexctmp = [Hexc,zeros(1,L-Ntmp)];
Hexctmp = reshape(Hexctmp,L,length(Hexctmp)/L);
for ii = 2:numel(Hexctmp)/L
    tmp = Hexctmp(:,ii);
    if sum(diff(tmp)) == 2*(L-1)
        Hexctmp(randi(L),ii) = 0;
    end
end
Hexc = Hexctmp(:);
Hexc(Hexc == 0) = [];


M    = length(Hexc);  % Number of excited harmonics

%% Generate frequency domain multisine I
I = zeros(Npp,1);
amplitudes = 1./(Hexc.^decay);         % Amplitudes with possible decay over frequency                         
phases     = 2*pi*rand(size(Hexc));    % Random uniformly distributed phases between 0 and 2\pi
I(Hexc+1) = amplitudes.*exp(1j*phases);            

%% Generate time domain signal i 
i = 2*real(ifft(I));                   % Inverse discrete Fourier transform
i = irms*i/rms(i);                     % Scale signal to have desired RMS value

I = fft(i)/Npp;

%% Show multisine
figure
plot(t,i,'.');
xlabel('time [s]','Interpreter','latex','FontSize',12);
ylabel('$i(t)$ [A]','Interpreter','latex','FontSize',12);
xlim([0,Tp]);


figure
loglog(f(1:end/2),abs(I(1:end/2)),'.');
hold on
stem(f(Hexc+1),abs(I(Hexc+1)),'x','MarkerSize',10,'color',colors(1));
set(gca,'xscal','log'); set(gca,'yscal','log');
grid on
xlabel('frequency [Hz]','Interpreter','latex','FontSize',11);
ylabel('$\vert I(\omega)\vert$ [A]','Interpreter','latex','FontSize',11);
xlim([fmin,fs/2]); 


%% Color function
function color = colors(i)
    MatlabColors = [0, 0.4470, 0.7410;
        0.8500, 0.3250, 0.0980];
    color = MatlabColors(i,:);
end




