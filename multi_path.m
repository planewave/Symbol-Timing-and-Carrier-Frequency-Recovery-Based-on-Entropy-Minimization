% multi path channel performance
clear
close all
clc

% EbN0 = (0:5:40)-4.77;
% num = length(EbN0);
M = 2; 
% M=2, var = 10e-5;  
% M=4, var=4.7 4.4 4.3
% M=8, var=4.6 4.1 4
% M=16 var=4
% 
k = log2(M); 
% EsN0 = EbN0 + 10*log10(k);
% EbN0 = 5:10:40;

% EbN0 = EsN0-10*log10(k);
% EsN0 = EbN0+ 10*log10(k);

% EsN0 = 5:5:40;
% EsN0_linear = 10.^(EsN0/10);
rolloff = 0.5;
T = 1;
L0 = 100;% 
% xi = 1/12+rolloff^2*(1/4-2/pi^2);
% mcrb_timing = 1./(8*pi^2*xi*L0*EsN0_linear)*T^2;

span = 10;       % Filter span
sps = 40;        % Samples per symbol = T

num = 80;% 
num_trial = num*L0;

%%
% rolloff = 0.75;
rrcFilter = rcosdesign(rolloff,span,sps);

% txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor',rolloff, ...
%     'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',sps);
% 
% rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor',rolloff, ...
%     'FilterSpanInSymbols',span,'InputSamplesPerSymbol',sps, ...
%     'DecimationFactor',1); % no decimation
data = randi([0 1],k*num_trial,1);
% hMod = comm.QPSKModulator('BitInput',true);
% hMod = comm.PSKModulator('BitInput',true,'ModulationOrder',M,'PhaseOffset',pi/M);
hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true);
modData = step(hMod, data);
% txSig = step(txfilter,modData);
txSig = upfirdn(modData, rrcFilter, sps)*sqrt(sps);
%%
Itr = 500;
errRateMid = zeros(1,Itr);
errRateOm = zeros(1,Itr);
errRateEm = zeros(1,Itr);
bSave={};
for itr = 1:Itr
    b = [1 zeros(1,randi([50 100])) 1.4*(rand-0.5) zeros(1,randi([0 50])) rand-0.5];
    % b = [zeros(1,5) 1];
    % b = [1 zeros(1,40) -0.5 ];
    % b = 1;
    txSig2 = filter(b,1,txSig);
    noisySig = awgn(txSig2,-0,'measured');
    
    rxFilt = upfirdn(noisySig, rrcFilter,1,1);
    sig = rxFilt(sps*10+1:end-sps*9-1);
    
    rxDownMid = sig(1:sps:end);
    rxDataMid = pskdemod(rxDownMid, 2, pi);
    errRateMid(itr) = sum(rxDataMid~=data)/num_trial;
    
    
    tauOm = round(-angle(sum(abs(sig).^2.*exp(-1j*2*pi*(0:length(sig)-1)'/sps)))/2/pi*sps);
    if tauOm>sps/2
        tauOm = tauOm-sps;
    end
    sigOm = rxFilt((sps*10+1:end-sps*9-1)+tauOm);
    rxDownOm =  sigOm(1:sps:end);
    rxDataOm = pskdemod(rxDownOm, 2, pi);
    errRateOm(itr) = sum(rxDataOm~=data)/num_trial;
    % rxDownMid = upfirdn(noisySig, rrcFilter,1,sps);
    % rxDownMid = rxDownMid(span+1:end-span);
    
    % close all
    % plotEye(sigOm(1:4001),sps,sps/2+1);
    % eyediagram(sig(1:8000),1*sps,sps)
    % grid
    % rxSig = step(rxfilter,noisySig);
    % sig = rxSig(12*sps+1:10000);
    
    
    txL0_rs = reshape(sig,[sps,num_trial]);
    entp = zeros(1,sps/2);
    for m = 1:sps/2
        if m>sps/4
            mm = m+20;
        else
            mm = m;
        end
        entp(m) = simpEntpEng(txL0_rs(mm,:),0.5,1);
    end
    [~, tauEm] = min(entp);
    if tauEm>sps/4
        tauEm = tauEm-sps/2+1;
    end
    
    sigEm = rxFilt((sps*10+1:end-sps*9-1)+tauEm);
    rxDownEm =  sigEm(1:sps:end);
    rxDataEm = pskdemod(rxDownEm, 2, pi);
    errRateEm(itr) = sum(rxDataEm~=data)/num_trial;
    if errRateEm(itr)-errRateOm(itr) > 0.01
        bSave{length(bSave)+1} = b;
    end

    if mod(itr,10) == 0
%         formatSpec = 'Iteration %d \n Error rate: Mid: %2.3f %, O&M: %2.3f %, EM: %2.3f %';
%         fprintf(formatSpec,itr,mean(errRateMid)*100,mean(errRateOm)*100,mean(errRateEm)*100)
        itr
    end
    
end
mean(errRateMid)
mean(errRateOm)
mean(errRateEm)
%%
deBug = false;


if deBug
%     b = bSave{11};
    b = [1, zeros(1,sps*1.7), -0.6, zeros(1,sps*0.4), 0.4];

    txSig2 = filter(b,1,txSig)/10;
    noisySig = awgn(txSig2, 0,'measured'); % EbN0 = SNR + 10*log10(sps) = 15
    
    rxFilt = upfirdn(noisySig, rrcFilter,1,1);
    sig = rxFilt(sps*10+1:end-sps*9-1);
    txL0_rs = reshape(sig,[sps,num_trial]);
    entp = zeros(1,sps);
    for m = 1:sps
        entp(m) = simpEntpEng(txL0_rs(m,:),0.5,0.5);
    end
    entp = circshift(entp,[0 (sps+0)/2]);
    % figure
    % t = -0.5:0.02:0.5-0.02;
    t = linspace(-0.5, 0.5, sps);
    
    close
    fg = figure;
    fg.Position=[600 -100 600 700];
%     subplot(2,1,1)
    plotEye(sig(1:3501),sps,sps/2+1);
    % plotEye(rxSig(12*sps+sps/2+1:15000),sps);
    ylim([-2 2])
    ax = gca;
    ax.XTick = -0.5:0.25:0.5;
    myplotsetting
%     subplot(2,1,2)
    % eyeX = (-sps/2:sps/2)'/sps;
%     plot(t,entp)
%     ylim([0.9 1])
%     ax = gca;
%     ax.XTick = -0.5:0.25:0.5;
    
    % legend('Modified', 'Original')
    xlabel('Normalized symbol timing offset')
    ylabel('Eye diagram entropy')
    myplotsetting
    
    clc
    
    [~, tauEm] = min(entp);
    sigEm = rxFilt((sps*10+1:end-sps*9-1)+tauEm-sps/2); % -sps/2 shift it back ...
    rxDownEm =  sigEm(1:sps:end);
    rxDataEm = pskdemod(rxDownEm, 2, pi);
    errRateEm = sum(rxDataEm~=data)/num_trial
    
    
    tauOm = round(-angle(sum(abs(sig).^2.*exp(-1j*2*pi*(0:length(sig)-1)'/sps)))/2/pi*sps)
    if tauOm>sps/2
        tauOm = tauOm-sps;
    end
    sigOm = rxFilt((sps*10+1:end-sps*9-1)+tauOm);
    rxDownOm =  sigOm(1:sps:end);
    rxDataOm = pskdemod(rxDownOm, 2, pi);
    errRateOm = sum(rxDataOm~=data)/num_trial
end
