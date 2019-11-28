% b1 is taps of RRC filter
% b2 is taps of RC filter
% b2=conv(b1,b1,'same');

clc;clear;close all;
rolloff=0.25;   % rolloff factor
span=6;             % time domain: [-span/2,span/2]
sps=4;          % sample per symbol
b1 = rcosdesign(rolloff,span,sps,'sqrt');
b2 = rcosdesign(rolloff,span,sps,'normal');

figure
t_b=-span/2:1/sps:span/2;
stem(t_b,b1,'b')
hold on 
stem(t_b,b2,'r')
legend('RRC','RC')

b=conv(b1,b1,'same');

figure
t_b=-span/2:1/sps:span/2;
stem(t_b,b,'b')
hold on 
stem(t_b,b2,'r')
legend('RRC','RC')