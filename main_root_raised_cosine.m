clc;clear;close all;
%% GNUradio RRC filter
samples_per_symbol=2;
excess_bw = 0.35;
nfilts = 32;
ntaps = nfilts * 11 * samples_per_symbol;    % make nfilts filters of ntaps each
[taps,t,num,den]=root_raised_cosine(nfilts,nfilts,1.0,excess_bw,ntaps);
figure
h1=plot(t,taps,'r');
title('root raised cosine')
grid on
hold on

excess_bw = 0.01;
[taps_0,t_0,num,den]=root_raised_cosine(nfilts,nfilts,1.0,excess_bw,ntaps);
h2=plot(t_0,taps_0,'g');

excess_bw = 1;
[taps_1,t_1,num,den]=root_raised_cosine(nfilts,nfilts,1.0,excess_bw,ntaps);
h3=plot(t_1,taps_1,'b');
legend([h2,h1,h3],'alpha=0.01','alpha=0.35','alpha=1')

%frequency scope
taps_f=fft(taps,length(taps));
taps_0_f=fft(taps_0,length(taps_0));
taps_1_f=fft(taps_1,length(taps_1));

figure
h1=plot(t,abs(fftshift(taps_f)),'r');
hold on
h2=plot(t_0,abs(fftshift(taps_0_f)),'g');
h3=plot(t_1,abs(fftshift(taps_1_f)),'b');
legend([h2,h1,h3],'alpha=0.01','alpha=0.35','alpha=1')

%% my RRC filter
alpha = 0.5;
fs = 20;  %sampling frequency
sequence=[-fs:1/fs:fs];
Num = (4*alpha*sequence).*cos((1+alpha)*pi*sequence)+sin((1-alpha)*pi*sequence);
Den = pi*sequence.*(1-(4*alpha*sequence).^2);
gt_alpha_root_05 = Num./Den;
DenZero_05 = find(abs(Den)<eps);
for i=1:length(DenZero_05)
    if sequence(DenZero_05(i))==0
        gt_alpha_root_05(DenZero_05(i))=1-alpha+4*alpha/pi;
    elseif abs(sequence(DenZero_05(i)))==1/(4*alpha)
        t=sequence(DenZero_05(i));
        num_gt_alpha_root=4*alpha*(cos(t*pi*(1+alpha))-t*pi*(1+alpha)*sin(t*pi*(1+alpha)))+(1-alpha)*pi*cos((1-alpha)*pi*t);
        den_gt_alpha_root=-32*pi*t^2*alpha^2;
        gt_alpha_root_05(DenZero_05(i))=num_gt_alpha_root/den_gt_alpha_root;
%         gt_alpha_root_05(DenZero_05(i))=alpha*((1+2/pi)*sin(pi/(4*alpha))+(1-2/pi)*cos(pi/(4*alpha)));
    end
end
N = length(gt_alpha_root_05);
GF_alpha_root_05 = fft(gt_alpha_root_05,N);

excess_bw = 0.5;
[taps_05,t_05,num,den]=root_raised_cosine(nfilts,nfilts,1.0,excess_bw,ntaps);

figure
plot(t_05,taps_05,'b','LineWidth',2);
hold on
plot(sequence,gt_alpha_root_05,'r','LineWidth',1)
grid on



