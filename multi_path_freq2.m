% multi path channel CFO
% borrow code from 'MCRB_freq_revisit2.m'
clear
close all
clc
gpu = false;
M = 4;
k = log2(M);
rolloff = 0.5;
T = 1;
L0 = 50;
span = 10;       % Filter span
sps = 40;        % Samples per symbol = T
num = 60;%
num_trial = num*L0;


%%
rrcFilter = rcosdesign(rolloff,span,sps);
data = randi([0 3],num_trial,1);
modData = pskmod(data,M,pi/M);
txSig = upfirdn(modData, rrcFilter, sps)*sqrt(sps);
%%
% Itr = 400;
% errRateMid = zeros(1,Itr);
% errRateOm = zeros(1,Itr);
% errRateEm = zeros(1,Itr);
% bSave={};
% for itr = 1:Itr
%     b = [1 zeros(1,randi([50 100])) ...
%         1.4*(rand-0.5) zeros(1,randi([0 50])) rand-0.5];
% %     b = [zeros(1,5) 1];
%     % b = [1 zeros(1,40) -0.5 ];
% %     b = 1;
%     txSig2 = filter(b,1,txSig)/10;
% %     txSig2 = circshift(txSig2,[-3, 0]);
%     noisySig = awgn(txSig2,-0,'measured');
%
%     rxFilt = upfirdn(noisySig, rrcFilter,1,1);
%     sig = rxFilt(sps*10+1:end-sps*9-1);
%
%     rxDownMid = sig(1:sps:end);
%     rxDataMid = pskdemod(rxDownMid, 2, pi);
%     errRateMid(itr) = sum(rxDataMid~=data)/num_trial;
%
%
%     tauOm = round(-angle(sum(abs(sig).^2.*exp ...
%         (-1j*2*pi*(0:length(sig)-1)'/sps)))/2/pi*sps);
%     if tauOm>sps/2
%         tauOm = tauOm-sps;
%     end
%     sigOm = rxFilt((sps*10+1:end-sps*9-1)+tauOm);
%     rxDownOm =  sigOm(1:sps:end);
%     rxDataOm = pskdemod(rxDownOm, 2, pi);
%     errRateOm(itr) = sum(rxDataOm~=data)/num_trial;
%     % rxDownMid = upfirdn(noisySig, rrcFilter,1,sps);
%     % rxDownMid = rxDownMid(span+1:end-span);
%
%     % close all
%     % plotEye(sigOm(1:4001),sps,sps/2+1);
%     % eyediagram(sig(1:8000),1*sps,sps)
%     % grid
%     % rxSig = step(rxfilter,noisySig);
%     % sig = rxSig(12*sps+1:10000);
%
%
%     txL0_rs = reshape(sig,[sps,num_trial]);
%     entp = zeros(1,sps/2);
%     for m = 1:sps/2
%         if m>sps/4
%             mm = m+20;
%         else
%             mm = m;
%         end
%         entp(m) = simpEntpEng(txL0_rs(mm,:),0.3,0.2);
%     end
%     [~, tauEm] = min(entp);
%     if tauEm>sps/4
%         tauEm = tauEm-sps/2-1; % '-1' for some reason
%     else
%         tauEm = tauEm-1;
%     end
%
%     sigEm = rxFilt((sps*10+1:end-sps*9-1)+tauEm);
%     rxDownEm =  sigEm(1:sps:end);
%     rxDataEm = pskdemod(rxDownEm, 2, pi);
%     errRateEm(itr) = sum(rxDataEm~=data)/num_trial;
%     if errRateEm(itr)-errRateOm(itr) > 0.01
%         bSave{length(bSave)+1} = b;
%     end
%
%     if mod(itr,10) == 0
%         itr
%         length(bSave)
%     end
%
% end
% mean(errRateMid)
% mean(errRateOm)
% mean(errRateEm)
%
%%
lenSig = length(txSig);
%     b = bSave{12};
% b = [1, zeros(1,sps*1.6), 0.6, zeros(1,sps*0.5), -0.3];
% b = 1;
% txSig2 = filter(b,1,txSig)/10; % multipath channel

carrier1 = exp(1j*(2*pi* rand * 0.005/sps *(1:lenSig)'));

txSig2 = 0.6*circshift(txSig,[sps*1.8, 0]);
carrier2 = exp(1j*(2*pi* rand*(0.01)/sps *(1:lenSig)'+3*rand-1));

txSig3 = 0.3*circshift(txSig,[sps*2.5, 0]);
carrier3 = exp(1j*(2*pi* rand*(0.01)/sps *(1:lenSig)'+1*rand-0.5));

dpMpSig = (txSig.*carrier1+txSig2.*carrier2+txSig3.*carrier3)/10;
% dpMpSig = txSig.*carrier1+txSig2.*carrier2;

noisySig = awgn(dpMpSig, 3,'measured'); 
% EbN0 = SNR + 10*log10(sps) = 15


% carrier = exp(1j*2*pi* 0.01/sps *(1:length(noisySig))');
% noisySig = carrier.*noisySig;


rxFilt = upfirdn(noisySig, rrcFilter,1,1); % matched filter
% sig = rxFilt(sps*10+1:end-sps*9-1); % truncation
% plotEye(sig(sps*600+1:sps*700+1),sps,sps/2+1);
%%
sigMid = rxFilt((sps*10+1:end-sps*9-1));

fs = 10;
% sigDown = sigMid;

[fEm, Data_comp, ~]= freqEstEntp(sigMid(1:sps:end) ...
    ,fs,0.4, 0.0, gpu);
rxFrqEm = Data_comp*exp(1j*0);
% if gpu 
%     rxFrqEm=gather(rxFrqEm);
% end
scatterplot(rxFrqEm)
% rxDataEm = pskdemod(rxFrqEm, M, pi/M);
% errRateEm = sum(rxDataEm~=data)/num_trial
fEm


[fML, Data_compML, ~]= freqEstML(sigMid(1:sps:end),fs,M);
rxFrqML = Data_compML*exp(-0.0j);
scatterplot(rxFrqML)
% rxDataML = pskdemod(rxFrqML, M, pi/M);
% errRateML = sum(rxDataML~=data)/num_trial
fML