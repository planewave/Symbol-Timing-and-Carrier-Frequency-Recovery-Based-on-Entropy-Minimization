% OM 前馈算法性能仿真 p405 图7.31 8psk 
%  
clear
close all
clc
L0 = 100;
EsN0 = 0:5:40;
EsN0_linear = 10.^(EsN0/10);
rolloff = 0.25; % 如果改成0.75 就是图7.29
T = 1;
% 计算 MCRB
xi = 1/12+rolloff^2*(1/4-2/pi^2);
mcrb_timing = 1./(8*pi^2*xi*L0*EsN0_linear)*T^2;
%%
span = 10;       % Filter span
sps = 4;        % Samples per symbol = T

% 这里没有修改默认rc幅度, 所以码元能量默认是1
rc = rcosdesign(rolloff, span, sps);
rc = rc/max(rc);
num = 600;% 仿真次数
L0 = 100;% 每次读取L0个码元
num_trial = num*L0;
%%
% modData = 2*randi([0 1],num_trial+100,1)-1; % BPSK
data = randi([0 1],3*num_trial+600,1);
% hMod = comm.QPSKModulator('BitInput',true);
hMod = comm.PSKModulator('BitInput',true,'ModulationOrder',8,'PhaseOffset',0);
modData = step(hMod, data);


txSig = upfirdn(modData, rc, sps); 
% txSig(1) = txSig(1)+1e-10j;
for p = 1:length(EsN0)

    rxSig = awgn(txSig,...
        EsN0(p)-10*log10(sps),'measured'); %
    rxSig = filter(rc,1,rxSig); % matched filter
    rxSig = rxSig(sps*20+1:end-sps*20);

    tau_OM = zeros(1,num);
    for n = 1:num
        txL0 = rxSig((n-1)*L0*sps+1:n*L0*sps);
  
        tau_OM(n) = angle(sum(txL0.*conj(txL0).*exp(-1j*2*pi*(0:sps*L0-1)'/sps)))/2/pi;
        
    end
    OM_var(p) = var(tau_OM);
end
%%
% semilogy(EsN0,mcrb_timing,EsN0,entp_var,'*-',EsN0,OM_var,'s-')
% legend('MCRB','EM','OM')
semilogy(EsN0,mcrb_timing,EsN0,OM_var,'s-','MarkerSize',12)
ax = gca;
ax.XTick = [0:5:40];
legend('MCRB','OM')
xlabel('E_s/N_0 (dB)')
ylabel('Timing varance')
ylim([1e-6 1e-1])
% ylim([1e-7 1e-2])
myplotsetting


% twit('跑完了 回来看! #matlab')