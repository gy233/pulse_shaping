% http://www.ni.com/tutorial/5657/en/
% https://wiki.gnuradio.org/index.php/Guided_Tutorial_PSK_Demodulation   CHAPTER: Adding a Channel
% Modulation Measurement Error Characterization
% For high order modulation schemes such as 64-QAM and 256-QAM
% small errors of DC offset, phase noise, quadrature skew, or IQ gain imbalance 
% can make the transitions of the RF signal too difficult to distinguish. 

clc;clear;close all;
Nsym=6;
beta=0.5;
sampsPerSym=64;
rctFilt = comm.RaisedCosineTransmitFilter(... 
    'Shape', 'Square root', ... 
    'RolloffFactor', beta, ... 
    'FilterSpanInSymbols', Nsym, ... 
    'OutputSamplesPerSymbol', sampsPerSym);
% b = coeffs(rctFilt);  % unit energy: sum(b.Numerator.^2)==1
% rctFilt.Gain = 1/max(b.Numerator); % Normalize to obtain maximum filter tap value of 1 


DataL = 100;                    % Data length in symbols, when DataL = 200, it means that there are 200 bits of data need to be sent
R = 1000;                       % Data rate, when R = 1000, it means that the system sends 1000 bit of data per second
Fs = R * sampsPerSym;           % Sampling frequency
Fr = 1*10^6;                    % carrier frequency
sampsPerSym_rf = 640;
Fs_rf = R * sampsPerSym_rf;  % sampling frequency of carrier

%% transmit part
%% real
% Generate random data 
x = 2*randi([0 1], DataL, 1)-1; 

% Time vector sampled at symbol rate in milliseconds 
tx = 1000 * (0: DataL - 1) / R;
% Filter 
yo = rctFilt([x; zeros(Nsym/2,1)]); 

% Filter group delay, since raised cosine filter is linear phase and 
% symmetric. 
fltDelay = Nsym / (2*R); 
% Correct for propagation delay by removing filter transients 
yo = yo(fltDelay*Fs+1:end); 
to = 1000 * (0: DataL*sampsPerSym - 1) / Fs; 
% Plot data. 
figure(1)
subplot(121)
stem(tx, x, 'kx'); hold on; 
% Plot filtered data. 
plot(to, yo, 'b-'); hold off; 
% Set axes and labels. 
axis([0 25 -1.7 1.7]); 
xlabel('Time (ms)'); 
ylabel('Amplitude'); 
legend('Transmitted Data', 'sqrt raised cosine Data', 'Location', 'southeast')
title('data (real part) after rrc filter')

%% imag
% Generate random data 
x1 = 2*randi([0 1], DataL, 1)-1; 

% Time vector sampled at symbol rate in milliseconds 
tx1 = 1000 * (0: DataL - 1) / R;

% Filter 
yo1 = rctFilt([x1; zeros(Nsym/2,1)]); 

% Filter group delay, since raised cosine filter is linear phase and 
% symmetric. 
fltDelay = Nsym / (2*R); 
% Correct for propagation delay by removing filter transients 
yo1 = yo1(fltDelay*Fs+1:end); 
to1 = 1000 * (0: DataL*sampsPerSym - 1) / Fs; 
% Plot data. 
figure(1)
subplot(122)
stem(tx1, x1, 'kx'); hold on; 
% Plot filtered data. 
plot(to1, yo1, 'b-'); hold off; 
% Set axes and labels. 
axis([0 25 -1.7 1.7]); 
xlabel('Time (ms)'); 
ylabel('Amplitude'); 
legend('Transmitted Data', 'sqrt raised cosine Data', 'Location', 'southeast')
title('data (imaginary part) after rrc filter')

%% frequency domain
figure(2)
N=length(yo);
% real
subplot(121)
yo_f_domain=fft(yo,N); 
A_yo_f_domain=abs(yo_f_domain);
f=(0:floor((N+1)/2)-1)*Fs/N;   
plot(f,A_yo_f_domain(1:floor((N+1)/2)));
xlabel('Frequency / Hertz');
ylabel('幅度');
axis([0 5000 -inf inf])
title('data (real part) after rrc filter')

% imag
subplot(122)
yo1_f_domain=fft(yo1,N); 
A_yo1_f_domain=abs(yo1_f_domain);
f=(0:floor((N+1)/2)-1)*Fs/N;   
plot(f,A_yo1_f_domain(1:floor((N+1)/2)));
xlabel('Frequency / Hertz');
ylabel('幅度');
axis([0 5000 -inf inf])
title('data (imaginary part) after rrc filter')

%%
figure(3)
scatter(yo,yo1,'b*')
legend('sample per symbol = 64')
title('constellation')

%% Direct (Homodyne) Upconversion
trf= 1000 * (0: DataL*sampsPerSym_rf - 1) / Fs_rf; 
yrf=interp1(to,yo,trf,'previous');
yrf1=interp1(to1,yo1,trf,'previous');
A_noise=0;
RF=yrf.*cos(2*pi*Fr*trf/1000)-yrf1.*sin(2*pi*Fr*trf/1000)+A_noise*randn(size(trf));
figure
plot(trf,RF)
axis([0 25 -1.7 1.7]); 
xlabel('Time (ms)'); 
ylabel('Amplitude'); 
title('RF')

%% receive part
rcrFilt = comm.RaisedCosineReceiveFilter(... 
    'Shape', 'Square root', ... 
    'RolloffFactor', beta, ... 
    'FilterSpanInSymbols', Nsym, ... 
    'InputSamplesPerSymbol', sampsPerSym_rf, ... 
    'DecimationFactor', 1); 
rcrFilt_verify = comm.RaisedCosineReceiveFilter(... 
    'Shape', 'Square root', ... 
    'RolloffFactor', beta, ... 
    'FilterSpanInSymbols', Nsym, ... 
    'InputSamplesPerSymbol', sampsPerSym, ... 
    'DecimationFactor', 1); 
%%
lo_phase_offset_re=0*pi/180;
lo_phase_offset_im=0*pi/180;
lo_frequency_offset_r=0;

%% real
yr=RF.*cos(2*pi*(Fr+lo_frequency_offset_r)*trf/1000+lo_phase_offset_re);
subplot(221)
plot(trf,yr)
xlabel('Time (ms)'); 
ylabel('Amplitude'); 
title('downconversion (real part)')

% this figure demonstrate why we use RRC filter for pulse shaping
% because RRC filter is a low-pass filter (LPF), Two RRC filter is a RC
% filter
% attention: the last Nsym/2 of RRC-filtered RF signal is NaN
yo_r = rcrFilt([yr'; zeros(Nsym*sampsPerSym_rf/2,1)]);
yo_r = yo_r(fltDelay*Fs_rf+1:end);
subplot(223)
plot(trf,yo_r)
hold on
% baseband signal (for varification purpose)
yo_r_verify = rcrFilt_verify([yo; zeros(Nsym*sampsPerSym/2,1)]);
yo_r_verify = yo_r_verify(fltDelay*Fs+1:end);
plot(to,yo_r_verify,'r')
% Down-sampled Carrier signal
yo_r_2 = rcrFilt_verify([yr(1:sampsPerSym_rf/sampsPerSym:end)'; zeros(Nsym*sampsPerSym/2,1)]);
yo_r_2 = yo_r_2(fltDelay*Fs+1:end);
plot(to,yo_r_2,'g')
stem(tx, x, 'kx');
xlabel('Time (ms)'); 
ylabel('Amplitude'); 
title('real part')
legend('Carrier signal after RRC filter','baseband signal after RRC filter','Down-sampled Carrier signal after RRC filter')

%% imag
yr1=-RF.*sin(2*pi*(Fr+lo_frequency_offset_r)*trf/1000+lo_phase_offset_im);
subplot(222)
plot(trf,yr1)
xlabel('Time (ms)'); 
ylabel('Amplitude'); 
title('downconversion (imaginary part)')

% this figure demonstrate why we use RRC filter for pulse shaping
% because RRC filter is a low-pass filter (LPF), when we convolve two RRC 
% filters together, we get a raised cosine filter
% attention: the last Nsym/2 of RRC-filtered RF signal is NaN
yo1_r = rcrFilt([yr1'; zeros(Nsym*sampsPerSym_rf/2,1)]);
yo1_r = yo1_r(fltDelay*Fs_rf+1:end);
subplot(224)
plot(trf,yo1_r)
hold on
% baseband signal (for varification purpose)
yo1_r_verify = rcrFilt_verify([yo1; zeros(Nsym*sampsPerSym/2,1)]);
yo1_r_verify = yo1_r_verify(fltDelay*Fs+1:end);
plot(to,yo1_r_verify,'r')
% Down-sampled Carrier signal
yo1_r_2 = rcrFilt_verify([yr1(1:sampsPerSym_rf/sampsPerSym:end)'; zeros(Nsym*sampsPerSym/2,1)]);
yo1_r_2 = yo1_r_2(fltDelay*Fs+1:end);
plot(to,yo1_r_2,'g')
stem(tx1, x1, 'kx');
xlabel('Time (ms)'); 
ylabel('Amplitude'); 
title('imaginary part')
legend('Carrier signal after RRC filter','baseband signal after RRC filter','Down-sampled Carrier signal after RRC filter')

%%
figure
scatter(yo_r,yo1_r,'b*')
title('constellation')

figure
scatter(yo_r_2,yo1_r_2,'b*')
title('constellation')