% for paper figure: Performance of three carrier frequency recovery algorithms
% 真大论文 频率回复性能
%  仿照书 p92 图 3.9 研究 Es/N0 关系
clear
close all
T = 1;
L0 = 100; % Number of input symbols
% EsN0 = -3:18;
EsN0 = 5:5:30;
EsN0_lin = 10.^(EsN0/10);
freq_var_mcrb = 3./(2*pi^2*L0^3*EsN0_lin);
gpu = false;% no gpu does not accelerate
rolloff = 0.5;
span = 10;
nSamp = 5;
txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',nSamp);

rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'InputSamplesPerSymbol',nSamp, ...
    'DecimationFactor',nSamp); % no decimation 

% xlim([-3 18])
% ylim([1e-8 1e-4])
%%
% sps = 4;           % Samples per symbol
% nSamp = nSym*sps;  % Number of samples
fs = 1;        % Sampling frequency (Hz)
% f_max = 1000;
f_off = 0.002;
% f_off = 0.01;

% pow = @(h) norm(h)^2/length(h);
freqOffset = comm.PhaseFrequencyOffset(...
    'FrequencyOffset',f_off, ...
    'SampleRate',fs);
% freqComp = comm.CoarseFrequencyCompensator(...
%     'Modulation','QPSK', ...
%     'SampleRate',fs, ...
%     'FrequencyResolution',1,...
%     'Algorithm','Correlation-based',...
%     'MaximumFrequencyOffset',f_max);
%%
num_trial = 500;
len_EsN0 = length(EsN0);
LnR_NDA = zeros(1,len_EsN0);
LnR_DA = zeros(1,len_EsN0);
Enp = zeros(1,len_EsN0);
Kay = zeros(1,len_EsN0);
gama = @(k) 3/2*L0/(L0^2-1)*(1-((2*k-L0)/L0).^2);
M=4;
r = 0.3:-0.05:0.05;
% r = 0.1*ones(1,6);
for m = 1:len_EsN0
    err_NDA = zeros(1,num_trial);
    err_DA = zeros(1,num_trial);
    err_Enp = zeros(1,num_trial);
    err_Kay = zeros(1,num_trial);
    parfor n = 1:num_trial
        data = randi([0 M-1],L0*2,1);
%         amp= randi([1 2],L0,1);
        modData = pskmod(data,M,pi/M);
%         txSig = modData;% 能量是1.5 
        txSig = step(txfilter,modData);
        SNR = EsN0(m) - 10*log10(nSamp);
        noisySig = awgn(txSig,SNR,'measured');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 看这里的噪声写法
%         noise = wgn(length(txSig),1,-EsN0(m),'complex');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rxSig = txSig+noise;
        rxSig = step(rxfilter,noisySig);
        rxSig = rxSig(51:150);
%         rxSig = awgn(txSig,EsN0(m),'measured');
        offsetData = step(freqOffset,rxSig);
        x = offsetData;
        order = M;
%         df = LuiseFreqEst(x,order,fs,f_max)
        raisedSignal = x .^ order;% 去调制
%         raisedSignal2 = x .* conj(modData);
    
    fs = 100;
    est_start = -1.1;
    entrp = zeros(1,21);
    eng = entrp;
    if gpu
        entrp = gpuArray(single(entrp));
    end
    for nn = 1:21
        f_tr = nn/10 + est_start;
        x_xp = x.*exp(-1j*2*pi*f_tr*(0:L0-1)/fs).';
        x_tr = raisedSignal.*exp(-1j*2*pi*f_tr*(0:L0-1)/fs).';
        entrp(nn) = simpEntpBat(x_xp,r(m),gpu);
        eng(nn) = abs(sum(x_tr));
    end
    [~,ind] = min(entrp);
    f_est = est_start+ind/10;
    
    [~,ind] = max(eng);
    f_est_ML = est_start+ind/10;
    
    for nn = 1:21
        f_tr = nn/100 + f_est-0.11;
        x_xp = x.*exp(-1j*2*pi*f_tr*(0:L0-1)/fs).';
        x_tr = raisedSignal.*exp(-1j*2*pi*f_tr*(0:L0-1)/fs).';
        entrp(nn) = simpEntpBat(x_xp,r(m),gpu);
        
        f_tr = nn/100 + f_est_ML-0.11;
        x_tr = raisedSignal.*exp(-1j*2*pi*f_tr*(0:L0-1)/fs).';
        eng(nn) = abs(sum(x_tr));
        
    end
    [~,ind] = min(entrp);
    f_est = f_est+ind/100-0.11;
    
    [~,ind] = max(eng);
    f_est_ML = f_est_ML+ind/100-0.11;
    
    
    for nn = 1:21
        f_tr = nn/1000 + f_est-0.011;
        x_xp = x.*exp(-1j*2*pi*f_tr*(0:L0-1)/fs).';
        x_tr = raisedSignal.*exp(-1j*2*pi*f_tr*(0:L0-1)/fs).';
        entrp(nn) = simpEntpBat(x_xp,r(m),gpu);
        
        f_tr = nn/1000 + f_est_ML-0.011;
        x_tr = raisedSignal.*exp(-1j*2*pi*f_tr*(0:L0-1)/fs).';
        eng(nn) = abs(sum(x_tr));
        
    end
    [~,ind] = min(entrp);
    f_est = f_est+ind/1000-0.011;    
    
    [~,ind] = max(eng);
    f_est_ML = f_est_ML+ind/1000-0.011;    
%     fs = 4;
%     f_est = freqEstEntp(offsetData, fs, 0.1, 0, gpu);
%     f_est_ML = freqEstML(offsetData, 1, M, 0);
    
    estFreqOP = 1/8/pi*angle(sum((x(2:end).*conj(x(1:end-1))).^4));
%         err_NDA(n) = estFreqOffset - f_off;
%         err_DA(n) = estFreqOffset2 - f_off;
        err_OP(n) = estFreqOP - f_off;
%         err_EM(n) = estFreqEM - f_off;
%         err_Kay(n) = v - f_off;
        err_EM(n) = f_est/fs- f_off;
        err_ML(n) = f_est_ML/M/fs- f_off;
    end
%     LnR_NDA(m) = var(err_NDA/fs); %Normalized frequency variance
%     LnR_DA(m) = var(err_DA/fs); 
    OP(m) = var(err_OP);
    EM(m) = var(err_EM);
    ML(m) = var(err_ML);
%     Kay(m) = var(err_Kay/fs);
%     Enp(m) = var(err_Enp/fs);
end
%%
close
fg = figure;
fg.Position=[600 200 600 500];
% semilogy(EsN0,freq_var_mcrb,EsN0,LnR_NDA,EsN0,LnR_DA,'-s',EsN0,Kay,'-+');
% legend('MCRB', 'L&R NDA','L&R DA','Kay')
% semilogy(EsN0,EM,'-.ok',EsN0,OP,'--dk',EsN0,ML,':^k',EsN0,freq_var_mcrb,'-k','MarkerSize',8);
semilogy(EsN0,EM,'-o',EsN0,OP,'-s',EsN0,ML,'-^',EsN0,freq_var_mcrb,'-','MarkerSize',12);
legend('EM','Open loop','ML','MCRB')


% ax = gca;
% ax.XTick = [-3:3:18];
% grid
xlabel('E_s/N_0 (dB)')
% xlim([-3 18])
ylim([1e-10 1e-4])
% grid
ylabel('Normalized frequency variance')
myplotsetting
grid

