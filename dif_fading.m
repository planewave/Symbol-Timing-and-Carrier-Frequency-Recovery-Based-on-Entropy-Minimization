% multi path channel performance
clear
close all
clc

gpu = true;

M = 4; 
k = log2(M); 
rolloff = 0.75;

span = 10;       % Filter span
sps = 40;        % Samples per symbol 
nSym = 3000; % don't explode the memory
baud = 240;
tSym = 1/baud;
train = 1000;

rrcFilter = rcosdesign(rolloff,span,sps);
data = randi([0 M-1],nSym,1);
hMod = comm.PSKModulator('ModulationOrder', M, 'PhaseOffset', pi/M);
hDmod = comm.PSKDemodulator('ModulationOrder', M, 'PhaseOffset', pi/M);
modData = step(hMod, data);
txSig = upfirdn(modData, rrcFilter, sps)*sqrt(sps);
%%
Itr = 500;
errRateMid = zeros(1,Itr);
errRateOm = zeros(1,Itr);
errRateEm = zeros(1,Itr);
chSave={};
chIf = zeros(1,Itr);
hStream   = RandStream.create('mt19937ar', 'seed', 12345);
parfor itr = 1:Itr

%     x = random('Exponential',1, 1, 2);
%     d1 = x(1);
%     d2 = x(1)+x(2);
%     a1 = random('Rayleigh', exp(-d1));
%     a2 = random('Rayleigh', exp(-d2));
%     ch = [a1, round(sps*d1), a2,round(sps*d2)];
%     chSave{itr} = ch;
%     txSig2 = ch(1)*circshift(txSig,[ch(2), 0]);
%     txSig3 = ch(3)*circshift(txSig,[ch(4), 0]);
%     mpSig = (txSig+txSig2+txSig3)/10;
%     noisySig = awgn(mpSig,3,'measured');

    chan = comm.RayleighChannel(...
        'SampleRate',baud*sps, ...
        'NormalizePathGains',true, ...
        'MaximumDopplerShift', 0, ...
        'PathGainsOutputPort',true);
    % RayleighChannel
    % RicianChannel
    
       chan.PathDelays = [0.0458 0.05264 0.05674 0.06357 0.07246 ]-0.0458;
       chan.AveragePathGains = 20*log10([7.737, 3.946, 4.848, 6.9 3.69 ]*1e-3);
%        chan.DirectPathDopplerShift = 3;

%         chan.PathDelays = [0 1.7*tSym];
%         chan.AveragePathGains = [0 -3];
%     chan.Visualization = 'Impulse response';
    [chanOut1,pathGains1] = chan(txSig);
%     plot(unwrap(angle(pathGains1)))

%     scatterplot(pathGains1(1:500:end, 1))
    
    noisySig = awgn(chanOut1,30);
%     scatterplot(noisySig)
    
%   matched filter
%   compensate for initial phase *conj(pathGains1(1,1))  if needed
    rxFilt = upfirdn(noisySig, rrcFilter,1,1); 
    sig = rxFilt(sps*span+1:end-sps*9-1); % compensate for the delay
    
    
    rxDownMid = sig(1:sps:end);
%     scatterplot(rxDownMid)
%     rxCar = carrSync(rxDownMid,modData);
%     rxEq = rxCar;
    
    % equalizer (one tap)
%     nWeights = 5;  % Single weight
    nFwdWeights = 6;  % Number of feedforward equalizer weights
    nFbkWeights = 6;  % Number of feedback filter weights
%     stepSize = 0.01; % Step size for LMS algorithm
    alg = rls(0.95);  % Adaptive algorithm object
    eqObj = dfe(nFwdWeights, nFbkWeights,alg,constellation(hMod).');
%     eqObj = lineareq(nWeights,alg,constellation(hMod).');  % Equalizer object
    [rxEq, rxEqd, err] = equalize(eqObj,rxDownMid, modData(1:train));
    
%     plot(abs(err))
%     figure
%     plot(rxEq(data==0),'o'); hold
%     plot(rxEq(data==1),'o');
%     plot(rxEq(data==2),'o');
%     plot(rxEq(data==3),'o');
%     xlim([-2 2]*1)
%     ylim([-2 2]*1)
%     hold off
    
%     scatterplot(rxEq)
    dmodMid = step(hDmod,rxEq);
    

    errRateMid(itr) = sum(dmodMid(train+1:end)~=data(train+1:end))/(nSym-train);
%     errRateMid(itr) = sum(dmodMid~=data)/nSym;
    
    tauOm = round(-angle(sum(abs(sig).^2.*exp(-1j*2*pi*(0:length(sig)-1)'/sps)))/2/pi*sps);
    if tauOm>sps/2
        tauOm = tauOm-sps;
    end
%     tauOm
    sigOm = rxFilt((sps*10+1:end-sps*9-1)+tauOm);
    rxDownOm =  sigOm(1:sps:end);
    eqObj.reset;
    [rxEq, rxEqd, err] = equalize(eqObj,rxDownOm, modData(1:train));
    rxDataOm = step(hDmod,rxEq);
%     errRateOm(itr) = sum(rxDataOm~=data)/nSym;
    errRateOm(itr) = sum(rxDataOm(train+1:end)~=data(train+1:end))/(nSym-train);
    % rxDownMid = upfirdn(noisySig, rrcFilter,1,sps);
    % rxDownMid = rxDownMid(span+1:end-span);
    
    % close all
    % plotEye(sigOm(1:4001),sps,sps/2+1);
    % eyediagram(sig(1:8000),1*sps,sps)
    % grid
    % rxSig = step(rxfilter,noisySig);
    % sig = rxSig(12*sps+1:10000);

%   EM 
    txL0_rs = circshift(reshape(sig,[sps,nSym]),[19, 0]);
%     txL0_rs = reshape(sig,[sps,numBits]);
    rxReshape = txL0_rs.';
    entp = bmre(rxReshape,0.20, 0.2, gpu);
%     entp = circshift(entp,[0 (sps+0)/2]);
    [~, tauEm] = min(entp(16:25));
    tauEm = tauEm-5;
%     tauEm
    

%     plot(rxReshape(:,20),'.')
    
    
%     entp = zeros(1,sps/2);
%     if gpu == true
%         entp = gpuArray(single(entp));
%     end
%     for m = 1:sps/2
%         if m>sps/4
%             mm = m+20;
%         else
%             mm = m;
%         end
%         entp(m) = simpEntpEng(txL0_rs(mm,:),0.3,0.2, gpu);
%     end
%     [~, tauEm] = min(entp);
%     if tauEm>sps/4
%         tauEm = tauEm-sps/2-1; % '-1' for some reason
%     else
%         tauEm = tauEm-1;
%     end
    
    sigEm = rxFilt((sps*10+1:end-sps*9-1)+tauEm);
    rxDownEm =  sigEm(1:sps:end);
    eqObj.reset;
    [rxEq, rxEqd, err] = equalize(eqObj,rxDownEm, modData(1:train));
    rxDataEm = step(hDmod,rxEq);
    errRateEm(itr) = sum(rxDataEm(train+1:end)~=data(train+1:end))/nSym;
%     if errRateEm(itr)<errRateOm(itr) 
%         errRateEm(itr)
%         errRateOm(itr)
%     end

    if mod(itr,25) == 0
%         formatSpec = 'Iteration %d \n Error rate: Mid: %2.3f %, O&M: %2.3f %, EM: %2.3f %';
%         fprintf(formatSpec,itr,mean(errRateMid)*100,mean(errRateOm)*100,mean(errRateEm)*100)
        itr
%         length(bSave)
    end
    
end
mid=mean(errRateMid)
om=mean(errRateOm)
em=mean(errRateEm)
(om-em)/om

chSave = chSave(chIf>0);

% plotEye(sig(sps*600+1:sps*700+1),sps,sps/2+1);
% hold
% 
% if gpu == true
%     tauEm = gather(tauEm);
% end
% 
% plot([1 1]*(tauEm)/sps, [-1.7 1.7],'Color',[0.8500    0.3250    0.0980])
% text((tauEm)/sps-0.03, -1.8, '\bf EM' )
% plot([1 1]*tauOm/sps, [-1.7 1.7],'Color',[0.8500    0.3250    0.0980])
% text((tauOm)/sps-0.04, -1.8, '\bf O&M' )
% myplotsetting
% hold off
%%

deBug = false;
% reset(gpuDevice(1))

if deBug
    close all
    
    ch = [0.5, round(sps*1.4), 0.2, round(sps*3.5)];
%     b1 = rand+1;
%     b2 = rand+2;
%     chn = [1, zeros(1,round(sps*1.4)), 0.5, zeros(1,round(sps*2.1)), 0.2];
%     chn = [1, zeros(1,round(sps*b1)), 0.5, zeros(1,round(sps*b2)), 0.2];
%     txSig2 = filter(chn,1,txSig)/10;
%     noisySig = awgn(txSig2, 0,'measured'); % EbN0 = SNR + 10*log10(sps) = 15
    
%     ch = chSave{2};
%     ch = [0.6*(rand), round(sps*(rand*3+1)), ...
%         0.3*(rand),round(sps*(rand*3+1))];

    txSig2 = ch(1)*circshift(txSig,[ch(2), 0]);
    txSig3 = ch(3)*circshift(txSig,[ch(4), 0]);
    mpSig = (txSig+txSig2+txSig3)/10;
    noisySig = awgn(mpSig, 0,'measured');


%     txSig4 = circshift(txSig,[02, 0]);
%     noisySig = awgn(txSig4/10, 0,'measured');
    
    rxFilt = upfirdn(noisySig, rrcFilter,1,1);
    sig = rxFilt(sps*10+1:end-sps*9-1);
    txL0_rs = circshift(reshape(sig,[sps,nSym]),[19, 0]);
    
    % fro mre fun
    rxReshape = txL0_rs.';
    entp = bmre(rxReshape,0.25, 0.1, gpu);
%     entp = circshift(entp,[0 (sps+0)/2]);
    [~, tauEm] = min(entp(16:25));
    tauEm = tauEm-5;
    
    
%     entp = zeros(1,sps);
%     for m = 1:sps
%         entp(m) = simpEntpEng(txL0_rs(m,:),0.3,0.2, gpu);
%     end
    
    
    
%     entp = circshift(entp,[0 (sps+0)/2]);
    % figure
    % t = -0.5:0.02:0.5-0.02;
%     t = linspace(-0.5, 0.5, sps);
    
%     close
%     fg = figure;
%     fg.Position=[600 -100 600 700];
%     subplot(2,1,1)
    plotEye(sig(sps*620+1:sps*700+1),sps,sps/2+1);
    % plotEye(rxSig(12*sps+sps/2+1:15000),sps);
    ylim([-2 2])
    ax = gca;
    ax.XTick = -0.5:0.25:0.5;
    
    hold
%     plot([0 0], [-1.7 1.7],'Color',[0.8500    0.3250    0.0980])
%     subplot(2,1,2)
    % eyeX = (-sps/2:sps/2)'/sps;
%     plot(t,entp)
%     ylim([0.9 1])
%     ax = gca;
%     ax.XTick = -0.5:0.25:0.5;
    
    % legend('Modified', 'Original')
%     xlabel('Normalized symbol timing offset')
%     ylabel('Eye diagram entropy')
%     myplotsetting
    
    clc
    
%     [~, tauEm] = min(entp);
%     tauEm = tauEm-sps/2; % don't why '-1', but has to be done
    if gpu
        tauEm = gather(tauEm);
    end
    sigEm = rxFilt((sps*10+1:end-sps*9-1)+tauEm); % -sps/2 shift it back ... 
    rxDownEm =  sigEm(1:sps:end);
    rxDataEm = pskdemod(rxDownEm, 2, pi);
    errRateEm = sum(rxDataEm~=data)/nSym
    plot([1 1]*(tauEm)/sps, [-1.7 1.7],'k')
%     plot([1 1]*(tauEm)/sps, [-1.7 1.7],'Color',[0.8500    0.3250    0.0980])
    text((tauEm)/sps-0.03, -1.8, '\bf EM' )    
    
    tauOm = round(-angle(sum(abs(sig).^2.*exp(-1j*2*pi*(0:length(sig)-1)'/sps)))/2/pi*sps);
    if tauOm>sps/2
        tauOm = tauOm-sps;
    end
    plot([1 1]*tauOm/sps, [-1.7 1.7],'k')
%     plot([1 1]*tauOm/sps, [-1.7 1.7],'Color',[0.8500    0.3250    0.0980])
    text((tauOm)/sps-0.04, -1.8, '\bf O&M' )
%     myplotsetting
    hold off
    
    sigOm = rxFilt((sps*10+1:end-sps*9-1)+tauOm);
    rxDownOm =  sigOm(1:sps:end);
    rxDataOm = pskdemod(rxDownOm, 2, pi);
    errRateOm = sum(rxDataOm~=data)/nSym
    
    
    sigMid = rxFilt((sps*10+1:end-sps*9-1));
    rxDownMid =  sigMid(1:sps:end);
    rxDataMid = pskdemod(rxDownMid, 2, pi);
    errRateMid = sum(rxDataMid~=data)/nSym
end
% figure
% plot(entp)
%% plot eye only
% close
% chn = [1, zeros(1,sps*1.6), 0.5, zeros(1,sps*0.4), 0.0];
% ch = chSave{1};
% txSig2 = ch(1)*circshift(txSig,[ch(2), 0]);
% txSig3 = ch(3)*circshift(txSig,[ch(4), 0]);
% mpSig = (txSig+txSig2+txSig3)/10;

% txSig2 = filter(chn,1,txSig)/10;
% noisySig = awgn(mpSig, 0,'measured');
%     
%     rxFilt = upfirdn(noisySig, rrcFilter,1,1);
%     sig = rxFilt(sps*10+1:end-sps*9-1);
%     plotEye(sig(sps*620+1:sps*720+1),sps,sps/2+1);