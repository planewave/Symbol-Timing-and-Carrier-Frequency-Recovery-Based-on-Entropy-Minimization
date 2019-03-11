%% 大论文 载波频偏 vs phase entropy 
clc
clear
close all
M = 4;                                   % Modulation order
k = log2(M);                             % Bits/symbol
n = 800;                                 % Transmitted bits
nSamp = 10;                              % Samples per symbol
EbNo = 15;                               % Eb/No (dB)

hMod = comm.QPSKModulator('BitInput',true);
hMod.PhaseOffset = pi/M;
span = 10;       % Filter span
rolloff = 0.5;
txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',nSamp);

rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'InputSamplesPerSymbol',nSamp, ...
    'DecimationFactor',nSamp); % no decimation 

data = randi([0 1],n,1);

modSig = step(hMod, data);
txSig = step(txfilter,modSig);

%% 载波频偏 
% edges = linspace(-pi,pi,32 );
f_res = 40; % 一共进行多少次 熵 估计
w_max = 2*pi/nSamp/100;
r = 0.3;
% f_offset_max = 0;
for n =1:f_res+1
    f_offset =(n-f_res/2-1)*w_max/f_res*2;
%     f_offset =(n-21)*0.0004/11;

    txSig_f = exp(-1j*f_offset*(1:length(txSig))').*txSig;
    
    SNR = EbNo + 10*log10(2) - 10*log10(nSamp);
    noisySig = awgn(txSig_f,SNR,'measured');
    rxSig = step(rxfilter,noisySig);
%     rcv = awgn(txSig3,15,'measured');
    
%     rcv = downsample(rcv,nSamp);
%     rcv_ang = angle(rcv);
%     [N,edges] = histcounts(rcv_ang,edges);
%     x = N/Ndata;
%     ent(n) = -sum(x.*log(eps+x));
    ent(n) = simpEntpBat(rxSig,r);
    ent2(n) = simpEntpBat(rxSig(1:50),r);
    ent3(n) = (simpEntpBat(rxSig(1:50),r)+ ...
        simpEntpBat(rxSig(50:100),r)+ ...
        simpEntpBat(rxSig(100:150),r)+ ...
        simpEntpBat(rxSig(150:200),r)+ ...
        simpEntpBat(rxSig(200:250),r)+ ...
        simpEntpBat(rxSig(250:300),r)+ ...
        simpEntpBat(rxSig(350:400),r))/7;
end
% histogram(rcv_ang,36)
%% plot here

tt = (-f_res/2:f_res/2)/(f_res/2);
plot(tt,ent,tt,ent2,tt,ent3)

% plot((-f_res/2:f_res/2)/(f_res/2),ent/max(ent))
xlabel('Normalized carrier frequency offset (%)')
% xlim([ -0.6 0.6])
% ax = gca;
% ax.XTick = [-0.6:0.2:0.6];
ylabel('Constellation diagram Entropy')
ylim([0.75 0.95])
legend('400 samples','50 samples','Block average')
grid
myplotsetting
%%
scatterplot(rxSig,1,10)
title('')
hold
scatter([0.707 0.707 -0.707 -0.707],[0.707 -0.707 0.707 -0.707],'r+')
ax = gca;
ax.XTick = [-1:0.5:1];
ax.YTick = [-1 -0.5 0 0.5 1];
axis equal
myplotsetting
hold off
close