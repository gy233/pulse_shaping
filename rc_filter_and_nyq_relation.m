% this program analysis Raised Cosine window, Nyquist impulse, RC filter
% and the realtionship between them. 
% book: Multirate Signal Processing for Communications Systems 
% chapter 4.1, 4.2 and 4.3
% surprising discovery: 
% Nyquist impulse spectrum power=sps/(2*t_range*f_symbol)
% I do not know why
clc;clear;close all;
f_symbol = 2;  % symbol frequency (number of symbols per second)
t_range = 10;
sps = 6;  % sample per symbol
fs = f_symbol*sps; % sample frequeny
t = -t_range:1/fs:t_range;

%% Raised Cosine window
alpha=[0.3 0.5 0.8];
% time domain
h_rcw_rec=[];
for i=1:length(alpha)
    h_rcw_num = cos(pi*alpha(i)*f_symbol*t);
    h_rcw_den = 1-(2*alpha(i)*f_symbol*t).^2;
    den_zero = find(abs(h_rcw_den)<=2*eps); % when f_symbol = 2; t_range = 2; sps = 6; abs(h_rcw_den) will equal to 2*eps
    h_rcw = h_rcw_num./h_rcw_den;
    h_rcw(den_zero) = pi*sin(pi*alpha(i)*f_symbol*t(den_zero))./(8*alpha(i)*f_symbol*t(den_zero));
    plot(t,h_rcw)
    hold on
    h_rcw_rec=[h_rcw_rec;h_rcw];
end
xlabel('Time (s)');
ylabel('Amplitude');
grid on
title('Raised Cosine window')
legend(strcat('alpha=',num2str(alpha(1))),strcat('alpha=',num2str(alpha(2))),strcat('alpha=',num2str(alpha(3))));

% frequency domain
figure(2)
for i=1:length(alpha)
    h_rcw=h_rcw_rec(i,:);
    N = length(h_rcw);
    Fs=fft(h_rcw,N);
    AFs=abs(fftshift(Fs)).^2/N;
    f=(0:N-1)*fs/N;
    plot(t*fs*fs/length(t),AFs);
    hold on
end
xlabel('Frequency (Hz)');
ylabel('Amplitude');
axis([-f_symbol,f_symbol,-inf,inf])
title('Raised Cosine window')
legend(strcat('alpha=',num2str(alpha(1))),strcat('alpha=',num2str(alpha(2))),strcat('alpha=',num2str(alpha(3))));

%% Nyquist impulse
% time domain
figure
h_nyq = sin(pi*f_symbol*t)./(pi*f_symbol*t);
h_nyq(t==0)=1;
plot(t,h_nyq)
xlabel('Time (s)');
ylabel('Amplitude');
grid on
title('Nyquist impulse')

% frequency domain
figure
N = length(h_nyq);
Fs=fft(h_nyq,N);
AFs=abs(fftshift(Fs)).^2/N;
f=(0:N-1)*fs/N;
plot(t*fs*fs/length(t),AFs);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Nyquist impulse')

%% RC filter
% time domain
figure
h_rcf_rec=[];
for i=1:length(alpha)
    h_rcw=h_rcw_rec(i,:);
    h_rcf = h_nyq.*h_rcw;
    plot(t,h_rcf)
    hold on
    h_rcf_rec=[h_rcf_rec;h_rcf];
end
xlabel('Time (s)');
ylabel('Amplitude');
grid on
title('RC filter')
legend(strcat('alpha=',num2str(alpha(1))),strcat('alpha=',num2str(alpha(2))),strcat('alpha=',num2str(alpha(3))));

% frequency domain
figure
for i=1:length(alpha)
    h_rcf=h_rcf_rec(i,:);
    N = length(h_rcf);
    Fs=fft(h_rcf,N);
%     AFs=abs(fftshift(Fs)).^2/sps.^2; % spectrum power density=1
    AFs=abs(fftshift(Fs)).^2/N;
    f=(0:N-1)*fs/N;
    plot(t*fs*fs/length(t),AFs);
    hold on
end
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('RC filter')
legend(strcat('alpha=',num2str(alpha(1))),strcat('alpha=',num2str(alpha(2))),strcat('alpha=',num2str(alpha(3))));
