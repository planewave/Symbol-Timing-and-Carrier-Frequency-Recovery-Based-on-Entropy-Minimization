% Journal paper Entropy Minimization Based 
% Symbol Timing and Carrier Frequency Recovery
%% Fig. 1
% plot eye diagram  
% Entropy Minimization Based Symbol Timing ...
close all
clear;
% nSym = 3e3;
% M = 2; 
% k = log2(M);                             % Bits/symbol
rolloff = 0.5;                           % Rolloff factor
% fc = 10e3;                               % carrier freq with offset
span = 6;                                % Filter span in symbols
% fsym = 1e3;                              % symbol rate
sps = 32;
% fs = sps*fsym;

NUM_WAV = 1;
% data = randi([0 M-1],nSym,1);
data2 = reshape(ones(500,1)*[0 1],[],1);
data = zeros(size(data2));
nSym = length(data);

loco = randperm(nSym);
for n = 1:nSym
    data(n) = data2(loco(n));
end


% modData = pskmod(data, M, pi/4);
modData = data*2-1;
% EbN0 = 0:10;                               % Eb/No (dB)
% EsN0 = EbN0 + 10*log10(k);                 % in case needed
% pulse shapingccc
rrcFilter = rcosdesign(rolloff,span,sps,'normal');
txPulse = upfirdn(modData, rrcFilter, sps)*sqrt(sps);% power = 1
% eyediagram(txPulse(sps*3+1:2000), sps)
sig = plotEye(txPulse(sps*3:5001), sps, sps/2+2);
hold
% plot([0 0], [-1.7 1.7],'Color',[0.8500    0.3250    0.0980])
plot([0 0], [-1.7 1.7],'--k')

text(0-0.02, -1.8, '\bf 1' )
% plot([0.125 0.125], [-1.7 1.7],'Color',[0.9290    0.6940    0.1250])
% plot([0.125 0.125], [-1.7 1.7],'Color',[0.8500    0.3250    0.0980])
plot([0.125 0.125], [-1.7 1.7],'--k')


text(0.125-0.02, -1.8, '\bf 2' )
% plot([0.25 0.25], [-1.7 1.7],'Color',[0.4940    0.1840    0.5560])
% plot([0.25 0.25], [-1.7 1.7],'Color',[0.8500    0.3250    0.0980])
plot([0.25 0.25], [-1.7 1.7],'--k')


text(0.25-0.02, -1.8, '\bf 3' )
% plot([0.375 0.375], [-1.7 1.7],'Color',[0.4660    0.6740    0.1880])
% plot([0.375 0.375], [-1.7 1.7],'Color',[0.8500    0.3250    0.0980])
plot([0.375 0.375], [-1.7 1.7],'--k')


text(0.375-0.02, -1.8, '\bf 4' )
% myplotsetting
hold off
% we will need histogram
nbins = 15;
ylimit = 0.6;
figure
subplot(2,2,1)
h=histogram(sig(17,:),nbins,'Normalization','probability');
% h.FaceColor = [0.500    0.50    0.5];
% h.EdgeColor = 'k';
ylim([0 ylimit])
ylabel('Probability')
xlabel({'Amplitude'})
title('1')
% myplotsetting

subplot(2,2,2)
h=histogram(sig(21,:),nbins,'Normalization','probability');
% h.FaceColor = [0.500    0.50    0.5];
% h.EdgeColor = 'k';
ylim([0 ylimit])
ylabel('Probability')
xlabel({'Amplitude'})
title('2')
% myplotsetting

subplot(2,2,3)
h=histogram(sig(25,:),nbins,'Normalization','probability');
% h.FaceColor = [0.500    0.50    0.5];
% h.EdgeColor = 'k';
ylim([0 ylimit])
ylabel('Probability')
xlabel({'Amplitude'})
title('3')
% myplotsetting

subplot(2,2,4)
h=histogram(sig(29,:),nbins,'Normalization','probability');
% h.FaceColor = [0.500    0.50    0.5];
% h.EdgeColor = 'k';
ylim([0 ylimit])
ylabel('Probability')
xlabel({'Amplitude'})
title('4')
% myplotsetting
%% Fig. 2
% constellation entropy demostration
clear
close all
% nSym = 50;
% M = 4;
nbin = 14;

% data = randi([0 M-1],nSym,1);
data2 = reshape(ones(14,1)*[0 1 2 3],[],1);
data = zeros(size(data2));
nSym = length(data);

loco = randperm(nSym);
for n = 1:nSym
    data(n) = data2(loco(n));
end
    


modData = pskmod(data, 4, pi/4);
subplot(3,2,1);
plot(modData,'.','MarkerSize',15)
grid
% xlabel({'Inphase','(a)-1'})
xlabel('Inphase')
ylabel('Quadrature')
title('(a)-1')
% scatterplot(modData, [],0,'.',h)
subplot(3,2,2)
h=histogram2(real(modData),imag(modData),nbin, ...
    'Normalization','probability');
% h.FaceColor = [0.500    0.50    0.5];
% h.EdgeColor = 'k';
xlim([-1 1])
ylim([-1 1])
zlim([0 0.26])
xlabel('Inphase')
ylabel('Quadrature')
zlabel('Probability')
title('(a)-2')
% axis equal
axesLabelsAlign3D
% 

modData2 = modData.*exp(1j*0.01*(1:nSym)');
subplot(3,2,3)
plot(modData2,'.','MarkerSize',15)
grid
% xlabel({'Inphase','(b)-1'})
xlabel('Inphase')
ylabel('Quadrature')
title('(b)-1')
subplot(3,2,4)
h=histogram2(real(modData2),imag(modData2),nbin, ...
    'Normalization','probability');
% h.FaceColor = [0.500    0.50    0.5];
% h.EdgeColor = 'k';
xlim([-1 1])
ylim([-1 1])
zlim([0 0.26])
xlabel('Inphase')
ylabel('Quadrature')
zlabel('Probability')
title('(b)-2')
axesLabelsAlign3D
% 

modData3 = modData.*exp(1j*0.02*(1:nSym)');
subplot(3,2,5)
plot(modData3,'.','MarkerSize',15)
grid
% xlabel({'Inphase','(c)-1'})
xlabel('Inphase')
ylabel('Quadrature')
title('(c)-1')
subplot(3,2,6)
h=histogram2(real(modData3),imag(modData3),nbin, ...
    'Normalization','probability');
% h.FaceColor = [0.500    0.50    0.5];
% h.EdgeColor = 'k';
xlim([-1 1])
ylim([-1 1])
zlim([0 0.26])
xlabel('Inphase')
ylabel('Quadrature')
zlabel('Probability')
title('(c)-2')
axesLabelsAlign3D
% scatterplot(modData2)



%% Not used
% ad-hoc entropy estimation
clear
close all
clc
L0 = 200;
EsN0 = -3:3:31;
num = length(EsN0);
EsN0_linear = 10.^(EsN0/10);


 % Create binary data for 4800, 2-bit symbols
       data = randi([0 1],2000,1);
 % Create a QPSK modulator System object with bits 
 %as inputs and Gray-coded signal constellation
       hMod = comm.QPSKModulator('BitInput',true);
 % Change the phase offset to pi/16
       hMod.PhaseOffset = pi/4;
 % Modulate and plot the data
       modData = step(hMod, data);
%
entp_histo= zeros(1,num);
entp_me = entp_histo;

r1 = 0.2;
r2 = 0.4;
for n = 1:num
    noiData = awgn(modData,EsN0(n));
    Xedges = -3:0.2:3;
    Yedges = Xedges;
    h = histogram2(real(noiData),imag(noiData),Xedges,Yedges);

    pn = h.Values/1000+eps;
    entp_histo(n) = -sum(sum((pn.*log2(pn))));
    
    Xedges2 = -3:0.4:3;
    Yedges2 = Xedges2;
    h2 = histogram2(real(noiData),imag(noiData),Xedges2,Yedges2);
    pn2 = h2.Values/1000+eps;
    entp_histo2(n) = -sum(sum((pn2.*log2(pn2))));
    
    entp_me(n) = simpEntpBat(noiData,r1);
    entp_me2(n) = simpEntpBat(noiData,r2);
    
end
%
close
fig=figure;
set(fig,'defaultAxesColorOrder',([[0 0 0]; [0 0 0]]));
yyaxis left
plot(EsN0,entp_me,'-ok',EsN0,entp_me2,'-dk','MarkerSize',8)
ylim([0.7 1.1])
ylabel('Modified R閚yi Entropy')

yyaxis right
plot(EsN0,entp_histo,'--^k',EsN0,entp_histo2,'--vk','MarkerSize',8)
ylim([1 9])
ylabel('Shannon Entropy (bit)')
xlabel('E_s/N_0 (dB)')
xlabel('SNR (dB)')
% myplotsetting

% myplotsetting
grid
xlim([-3 30])

legend('Modified R閚yi Entropy \it{r}=0.2','Modified R閚yi Entropy \it{r}=0.4', ...
    'Shannon Entropy bin=0.2','Shannon Entropy bin=0.4')


% entp_histoN = entp_histo/25.6+0.672;
% entp_histo2N = entp_histo2/25.6+0.672;
% plot(EsN0,entp_histoN,'-dk',EsN0,entp_histo2N,'-ok', ...
%     EsN0,entp_me,'--vk',EsN0,entp_me2,'--^k', ...
%     'MarkerSize',8);
% legend('Histogram bin=0.2','Histogram bin=0.4', ...
%     'MRE \it{r_{ag}}=0.2','MRE \it{r_{ag}}=0.4')
% ylim([0.7 1.05])
% xlabel('SNR (dB)')
% ylabel('Normalized Entropy')
% myplotsetting
% grid
% xlim([-3 30])


%% Fig. 3 and Fig. 4
% eye diagram entropy plot
clear
clc
close all
figure

M = 4;                                   % Modulation order
k = log2(M);                             % Bits/symbol
n = 1000;                                % Transmitted bits
nSamp = 40;                              % Samples per symbol
EbNo = 15;                               % Eb/No (dB)

hMod = comm.QPSKModulator('BitInput',true);
hMod.PhaseOffset = pi/M;
span = 10;       % Filter span
rolloff = 0.5;
txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',nSamp);

rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'InputSamplesPerSymbol',nSamp, ...
    'DecimationFactor',1); % no decimation 

data = randi([0 1],n,1);

modSig = step(hMod, data);
txSig = step(txfilter,modSig);

t = 0:1/40:500-1/40;
% uncomment below for Fig. 4
% txSig = txSig.*exp(1j*0.01*pi*t');

% eyediagram(txSig(1+nSamp*50:3000),nSamp)
SNR = EbNo + 10*log10(k) - 10*log10(nSamp);
noisySig = awgn(txSig,SNR,'measured');
% eyediagram(noisySig(1+nSamp*50:3000),nSamp)
rxSig = step(rxfilter,noisySig);
% eyediagram(rxSig(1+nSamp*50:3000),nSamp)

%
sps = nSamp;

eyeX = (-sps/2:sps/2)'/sps;
eyeData = real(rxSig(1+sps*25:sps*350));
nRep = length(eyeData)/sps; % 必须是整数!


% eyediagram(real(txSig(1+sps*25:sps*350)),sps)

sig_rs = reshape(rxSig(1+16.5*sps:416.5*sps),[sps,400]);
sig_rs = [sig_rs;sig_rs(1,:)];
% plot(real(sig_rs))


r = 0.3;
for n = 1:sps+1
    entp(n) = simpEntpBat(sig_rs(n,:),r);
    entp2(n) = simpEntpEng(sig_rs(n,:),r,0.3);
end

close
fg = figure;
fg.Position=[600 200 600 700];
subplot(2,1,1)

% yyaxis left

plotEye(rxSig(1+sps*55.5:sps*120),sps);
ylim([-1.45 1.5])
ax = gca;
ax.XTick = -0.5:0.25:0.5;
grid on
% set(gca,'xtick',[])
xlabel([])

myplotsetting
subplot(2,1,2)
% yyaxis right

plot(eyeX,entp2,'--',eyeX,entp,'-')
ylim([0.7 .99])
ax = gca;
ax.XTick = -0.5:0.25:0.5;

legend('BMRE', 'MRE')
xlabel('Normalized symbol period')
ylabel('Eye diagram entropy')
grid on
myplotsetting



