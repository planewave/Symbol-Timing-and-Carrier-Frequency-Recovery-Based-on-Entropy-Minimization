% for paper NDA  Entropy QPSK rolloff

clear
close all
clc
M = 4; 

k = log2(M); 
% EsN0 = EbN0 + 10*log10(k);
EbN0 = 5:5:40;

% EbN0 = EsN0-10*log10(k);
% EsN0 = EbN0+ 10*log10(k);

EsN0 = 5:5:40;
EsN0_linear = 10.^(EsN0/10);
rolloff = 0.25;
rolloff2 = 0.05;
T = 1;
L0 = 100;
xi = 1/12+rolloff^2*(1/4-2/pi^2);
mcrb_timing = 1./(8*pi^2*xi*L0*EsN0_linear)*T^2;

% xi2 = 1/12+rolloff2^2*(1/4-2/pi^2);
% mcrb_timing2 = 1./(8*pi^2*xi2*L0*EsN0_linear)*T^2;
% 
% semilogy(EsN0, mcrb_timing, EsN0, mcrb_timing2)

span = 10;       % Filter span
sps = 8;        % Samples per symbol = T

num = 800;% 仿真次数
num_trial = num*L0;
rrcFilter = rcosdesign(rolloff,span,sps);
data = randi([0 M-1],num_trial+600,1);

txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',sps);

rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'InputSamplesPerSymbol',sps, ...
    'DecimationFactor',1); % no decimation

txfilter2 = comm.RaisedCosineTransmitFilter('RolloffFactor',rolloff2, ...
    'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',sps);

rxfilter2 = comm.RaisedCosineReceiveFilter('RolloffFactor',rolloff2, ...
    'FilterSpanInSymbols',span,'InputSamplesPerSymbol',sps, ...
    'DecimationFactor',1); % no decimation


% data = randi([0 1],k*num_trial+600,1);
% data2 = randi([0 1],2*num_trial+600,1);
% hMod2 = comm.QPSKModulator('BitInput',true);
% hMod = comm.PSKModulator('BitInput',true,'ModulationOrder',M,'PhaseOffset',pi/M);
% hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true);


hMod = comm.PSKModulator('ModulationOrder', M, 'PhaseOffset', pi/M);
hDmod = comm.PSKDemodulator('ModulationOrder', M, 'PhaseOffset', pi/M);
modData = step(hMod, data);
% modData2 = step(hMod2, data2);
txSig = step(txfilter,modData);
txSig2 = step(txfilter2,modData);
% 加入频偏
% t = 0:1/sps:length(txSig)/sps-1/sps;
% txSig = txSig.*exp(1j*0.00*pi*t');
%%
r = [0.5:-0.05:0.15]; 
r = [0.55:-0.05:0.2]; 
% r = 0.5:-0.04:0.01;
% r = 0.3*ones(1,9);
for p = 1:length(EsN0)
%     SNR = 40 - 10*log10(sps);
    SNR = EsN0(p) - 10*log10(sps);
%     SNR = 16-10*log10(sps);
%     SNR = 5;
    noisySig = awgn(txSig,SNR,'measured');
    rxSig = step(rxfilter,noisySig);
    
    noisySig2 = awgn(txSig2,SNR,'measured');
    rxSig2 = step(rxfilter2,noisySig2);

%     rxSig = awgn(txSig,EbN0(p)-10*log10(sps),'measured');
%     rxSig = filter(rc,1,rxSig); % matched filter
%     rxSig = rxSig(sps*20+1:end-sps*20);
    
%%
    
    tau_entrp = zeros(1,num);
    tau_OM = tau_entrp;
    offset = 0;
    for n = 1:num 
        txL0 = rxSig((n+0)*L0*sps+1+offset:(n+1)*L0*sps+offset);
        txL0_rs = reshape(txL0,[sps,L0]);
        
        txL02 = rxSig2((n+0)*L0*sps+1+offset:(n+1)*L0*sps+offset);
        txL0_rs2 = reshape(txL02,[sps,L0]);
        
        entp = zeros(sps,1);
        entp2 = zeros(sps,1);
        for m = 1:sps
            entp(m) = simpEntpEng(txL0_rs(m,:),r(p),0.3);
            entp2(m) = simpEntpEng(txL0_rs2(m,:),r(p),0.3);
        end
%         txL0_rs = reshape(txL0,[sps*2,L0/2]);
%         entp = zeros(sps*2,1);
%         for m = 1:sps*2
%             entp(m) = simpEntpEng(txL0_rs(m,:),0.25,0.25);
%         end
        
        entp = circshift(entp,[(sps+0)/2 0]);
        entp2 = circshift(entp2,[(sps+0)/2 0]);
        
%         f = fit([-0.3:0.1:0.3]',entp(3:9),'poly2');
%         tau_entrp(n) = -f.p2/f.p1/2;

        tau_entrp(n) = angle(sum(entp.*exp(-1j*2*pi*(0:sps-1)'/sps)))/2/pi;
        tau_entrp2(n) = angle(sum(entp2.*exp(-1j*2*pi*(0:sps-1)'/sps)))/2/pi;
%         tau_entrp(n) = angle(sum(entp(1:10).*exp(-1j*2*pi*(0:sps-1)'/sps)))/2/pi;
        tau_OM(n) = angle(sum(abs(txL0).^2.*exp(-1j*2*pi*(0:sps*L0-1)'/sps)))/2/pi;
        tau_OM2(n) = angle(sum(abs(txL02).^2.*exp(-1j*2*pi*(0:sps*L0-1)'/sps)))/2/pi;
%         tau_OM(n) = angle(sum(txL0.*conj(txL0).*exp(-1j*2*pi*(0:sps*L0-1)'/sps)))/2/pi;
      
    end
%     entp_men(p) = mean(tau_entrp);
    entp_var(p) = var(tau_entrp);
    OM_var(p) = var(tau_OM);
    entp_var2(p) = var(tau_entrp2);
    OM_var2(p) = var(tau_OM2);
%     if p ==3
%         pause();
%     end
end
%% 比较 bpsk qpsk
fg = figure;
fg.Position=[600 200 600 500];
% semilogy(EsN0,entp_var,'--ok', ...
%     EsN0,OM_var,'--dk', ...
%     EsN0,entp_var2,'-.ok', ...
%     EsN0,OM_var2,'-.dk' ...
%     ,EsN0,mcrb_timing,'-k','MarkerSize',8);
semilogy(EsN0,entp_var,'-o', ...
    EsN0,OM_var,'-s', ...
    EsN0,entp_var2,'-.o', ...
    EsN0,OM_var2,'-.s' ...
    ,EsN0,mcrb_timing,'-','MarkerSize',12);
legend('\alpha=0.25, EM','\alpha=0.25, O&M','\alpha=0.05, EM','\alpha=0.05, O&M','MCRB');

% leg1 = legend('$\alpha=0.25, EM$','$\alpha=0.25, O\&M$','$\alpha=0.05, EM$','$\alpha=0.05, O\&M$','MCRB');
% set(leg1,'Interpreter','latex');

xlabel('E_s/N_0 (dB)')
ylabel('Normalized timing variance')
% ylim([10e-7 1.2e-2])
% ylim([10e-7 1.2e-2])
xlim([5 40])
ylim([10e-7 10e-2])
% title(['M = ' num2str(M) ', \alpha = ' num2str(rolloff)])
myplotsetting
grid
%%
% fg=figure;
% fg.Position=[600 200 600 550];
% semilogy(EsN0,entp_var,'o-', ...
%     EsN0,OM_var,'s-' ...
%     ,EsN0,mcrb_timing,'MarkerSize',12);
% 
% legend('EM','O&M','MCRB')
% xlabel('E_s/N_0 (dB)')
% ylabel('Normalized timing varance')
% ylim([1e-6 1.1e-2])
% % title(['M = ' num2str(M) ', \alpha = ' num2str(rolloff)])
% myplotsetting
%% 
% plot UWA channel  
figure
load('sample_10km_cir.mat')
plot(tau_x, abs(cir_first_idx(:,1)/max(cir_first_idx(:,1))))
set(get(gca,'Children'),'linewidth',1.5); % plot line
set(gca, 'LineWidth', 1.5); % frame line
xlabel('Time (s)')
ylabel('Normalized amplitude')
set(gca,'FontSize',15);
xlim([0 0.3])
% myplotsetting