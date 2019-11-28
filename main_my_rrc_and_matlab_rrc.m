clc;clear;close all;

% matlab tool
rolloff=0.35;   % rolloff factor
span=11*2;             % time domain: [-span/2,span/2]
sps=32;          % sample per symbol
shape='sqrt';   % 'sqrt' for root raised cosine or 'normal' for raised cosine
b = rcosdesign(rolloff,span,sps,shape);
b1 = rcosdesign(0.5,span,sps,shape);
b2 = rcosdesign(rolloff,span,sps,'normal');
b3 = rcosdesign(0.5,span,sps,'normal');

figure
t_b=-span/2:1/sps:span/2;
plot(t_b,b,'b')
hold on 
plot(t_b,b1,'r')
plot(t_b,b2,'g--')
plot(t_b,b3,'y--')
grid on
title('matlab rcosdesign')
legend('rolloff=0.35, shape=sqrt','rolloff=0.5, shape=sqrt','rolloff=0.35, shape=normal','rolloff=0.5, shape=normal')

% you can see span*sps samples in the impulse response via fvtool
% figure is the same as stem(b)
% fvtool(b,'Analysis','impulse');
% figure
% stem(t_b,b)
% title('taps of RRC filter, rolloff=0.35')

% pulse shaping
d = [1 1 1 1 1 -1 1 -1];
x = upfirdn(d, b, sps); % process of convolution
delay=(span/2)*sps;

% figure is the same as stem(x)
% fvtool(x,'Analysis','impulse');
figure
t_x = 0:1/sps:length(d)+span/2-1;
stem(t_x,x(1+delay:end))


%% my root raised cosine (transmit part)
symbol_rate=1;
sampling_freq=sps*symbol_rate;
gain=sampling_freq;  % in GNU Radio, gain==sampling_freq, witch make max(taps)==1; 
alpha=rolloff;
ntaps=span*sps;
[taps,t_b,num,den] = root_raised_cosine(gain,sampling_freq,symbol_rate,alpha,ntaps);
taps=taps/sqrt(sum(taps.^2)); % power of taps: sum(taps.^2)==1 

figure
stem(t_b,taps,'b','LineWidth',2)
hold on 
stem(t_b,b,'r','LineWidth',1)
legend('my rrc','matlab rrc (rcosdesign)')
title('taps of RRC filter')

% convolution
x_my_rrc=zeros(1,sps*(length(d)-1)+bitor(ntaps,1));
for i=1:length(d)
    x_my_rrc(sps*i-sps+1:sps*i-sps+bitor(ntaps,1))= x_my_rrc(sps*i-sps+1:sps*i-sps+bitor(ntaps,1))+taps*d(i);
end
x_my_rrc_delay=x_my_rrc(1+delay:end);
figure
plot(t_x,x_my_rrc_delay,'b','LineWidth',3)
hold on 
plot(t_x,x(1+delay:end),'r','LineWidth',2)
grid on
legend('my rrc filter','matlab filter (rcosdesign)')
title('data after RRC filter transmit')

%% my root raised cosine and comm.RaisedCosineTransmitFilter & comm.RaisedCosineReceiveFilter
x_my_rrc_r=conv(x_my_rrc_delay,taps,'same'); 

rctFilt = comm.RaisedCosineTransmitFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', rolloff, ...
    'FilterSpanInSymbols', span, ...
    'OutputSamplesPerSymbol', sps);

yo = rctFilt([d'; zeros(span/2,1)]);
yo=yo(1+delay:end);
t_yo=0:1/sps:length(d)-1/sps;
figure
plot(t_x,x_my_rrc_delay,'b','LineWidth',3)
hold on 
plot(t_yo,yo,'r')
legend('my rrc filter','matlab filter (comm.RaisedCosineTransmitFilter)')
title('data after RRC filter transmit')

rcrFilt = comm.RaisedCosineReceiveFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          rolloff, ...
  'FilterSpanInSymbols',    span, ...
  'InputSamplesPerSymbol',  sps, ...
  'DecimationFactor',       1);

yr = rcrFilt([yo; zeros(span*sps/2, 1)]);
yr=yr(1+delay:end);

figure
plot(t_x,x_my_rrc_r,'b','LineWidth',3)
hold on
plot(t_yo,yr,'r')
stem(0:length(d)-1,d)
legend('my rrc filter','matlab filter (comm.RaisedCosineReceiveFilter)','data')
title('conv(RRC,RRC,same) == RC')

%% if x_my_rrc=x_my_rrc(1+delay:end); doesn't delay, x_my_rrc_r will coincide with the first d
x_my_rrc_r=conv(x_my_rrc,taps,'same'); 
figure
plot(t_x,x_my_rrc_r(1+delay:end),'b','LineWidth',3)
hold on
plot(t_yo,yr,'r')
stem(0:length(d)-1,d)
legend('my rrc filter','matlab filter (comm.RaisedCosineReceiveFilter)','data')
title('conv(RRC,RRC,same) == RC')