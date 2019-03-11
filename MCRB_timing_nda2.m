% 性能仿真 NDA  Entropy 算法 
% 红皮书 p404 比较
clear
close all
clc

% EbN0 = (0:5:40)-4.77;
% num = length(EbN0);
M = 4; 
% M=2, var = 10e-5;  
% M=4, var=4.7 4.4 4.3
% M=8, var=4.6 4.1 4
% M=16 var=4
% 
k = log2(M); 
% EsN0 = EbN0 + 10*log10(k);
EbN0 = 5:10:40;

% EbN0 = EsN0-10*log10(k);
EsN0 = EbN0+ 10*log10(k);

EsN0 = 5:5:40;
EsN0_linear = 10.^(EsN0/10);
rolloff = 0.25;
T = 1;
L0 = 100;% 每次读取L0个码元
% 计算 MCRB
xi = 1/12+rolloff^2*(1/4-2/pi^2);
mcrb_timing = 1./(8*pi^2*xi*L0*EsN0_linear)*T^2;

span = 10;       % Filter span
sps = 10;        % Samples per symbol = T

num = 500;% 仿真次数
num_trial = num*L0;

%%
% rolloff = 0.75;
txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',sps);

rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'InputSamplesPerSymbol',sps, ...
    'DecimationFactor',1); % no decimation
data = randi([0 1],k*num_trial+600,1);
% hMod = comm.QPSKModulator('BitInput',true);
% hMod = comm.PSKModulator('BitInput',true,'ModulationOrder',M,'PhaseOffset',pi/M);
hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true);
modData = step(hMod, data);
txSig = step(txfilter,modData);
% 加入多径
b = [1 0 0 0 0 0.5-0.2j 0 0 -0.3];
txSig = filter(b,1,txSig);
% 加入频偏
% t = 0:1/sps:length(txSig)/sps-1/sps;
% txSig = txSig.*exp(1j*0.00*pi*t');
%%
r = [0.55:-0.05:0.15]; %熵距离参数
% r = 0.5:-0.04:0.01;
% r = 0.3*ones(1,9);
for p = 1:length(EsN0)
%     SNR = 40 - 10*log10(sps);
    SNR = EsN0(p) - 10*log10(sps);
%     SNR = 16-10*log10(sps);
%     SNR = 5;
    noisySig = awgn(txSig,SNR,'measured');
    rxSig = step(rxfilter,noisySig);

%     rxSig = awgn(txSig,EbN0(p)-10*log10(sps),'measured');
%     rxSig = filter(rc,1,rxSig); % matched filter
%     rxSig = rxSig(sps*20+1:end-sps*20);
    
    %%
    
    tau_entrp = zeros(1,num);
    tau_OM = tau_entrp;
    offset = 0;
    for n = 1:num % 500次测试
        txL0 = rxSig((n+0)*L0*sps+1+offset:(n+1)*L0*sps+offset);
        txL0_rs = reshape(txL0,[sps,L0]);
        entp = zeros(sps,1);
        for m = 1:sps
            entp(m) = simpEntpEng(txL0_rs(m,:),0.4,0.3);
        end
%         txL0_rs = reshape(txL0,[sps*2,L0/2]);
%         entp = zeros(sps*2,1);
%         for m = 1:sps*2
%             entp(m) = simpEntpEng(txL0_rs(m,:),0.25,0.25);
%         end
        
        entp = circshift(entp,[(sps+0)/2 0]);
        
%         f = fit([-0.3:0.1:0.3]',entp(3:9),'poly2');
%         tau_entrp(n) = -f.p2/f.p1/2;

        tau_entrp(n) = angle(sum(entp.*exp(-1j*2*pi*(0:sps-1)'/sps)))/2/pi;
%         tau_entrp(n) = angle(sum(entp(1:10).*exp(-1j*2*pi*(0:sps-1)'/sps)))/2/pi;
        tau_OM(n) = angle(sum(abs(txL0).^2.*exp(-1j*2*pi*(0:sps*L0-1)'/sps)))/2/pi;
%         tau_OM(n) = angle(sum(txL0.*conj(txL0).*exp(-1j*2*pi*(0:sps*L0-1)'/sps)))/2/pi;
      
    end
    entp_men(p) = mean(tau_entrp);
    entp_var(p) = var(tau_entrp);
    OM_var(p) = var(tau_OM);
%     if p ==3
%         pause();
%     end
end
eyediagram(rxSig(12*sps+1:1000),2*sps)
%% 比较 bpsk qpsk
%  load('timing_b.mat')
% semilogy(EbN0,entp_var,'o-', ...
%     EbN0,OM_var,'s-', ...
%     EbN0,entp_var_b,'o-.', ...
%     EbN0,OM_var_b,'s-.' ...
%     ,EbN0,mcrb_timing,'MarkerSize',12);
% legend('QPSK, EM','QPSK, O&M','BPSK, EM','BPSK, O&M','MCRB')
% xlabel('E_b/N_0 (dB)')
% ylabel('Normalized timing varance')
% ylim([1e-7 1e-1])
% % xlim([5 35])
% % title(['M = ' num2str(M) ', \alpha = ' num2str(rolloff)])
% myplotsetting
%%
fg=figure;
fg.Position=[600 200 600 550];
semilogy(EsN0,entp_var,'o-', ...
    EsN0,OM_var,'s-' ...
    ,EsN0,mcrb_timing,'MarkerSize',12);

legend('EM','O&M','MCRB')
xlabel('E_s/N_0 (dB)')
ylabel('Normalized timing variance')
ylim([1e-7 1e-1])
% title(['M = ' num2str(M) ', \alpha = ' num2str(rolloff)])
myplotsetting
