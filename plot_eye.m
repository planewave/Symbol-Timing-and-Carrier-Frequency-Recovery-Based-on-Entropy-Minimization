% plot eye diagram for paper
% Entropy Minimization Based Symbol Timing ...
% abolished
close all
clear;
nSym = 6e3;
M = 2; 
k = log2(M);                             % Bits/symbol
rolloff = 0.5;                           % Rolloff factor
fc = 10e3;                               % carrier freq with offset
span = 6;                                % Filter span in symbols
fsym = 1e3;                              % symbol rate
sps = 32;
fs = sps*fsym;

NUM_WAV = 1;
data = randi([0 M-1],nSym,1);
% modData = pskmod(data, M, pi/4);
modData = data*2-1;
EbN0 = 0:10;                               % Eb/No (dB)
EsN0 = EbN0 + 10*log10(k);                 % in case needed
% pulse shapingccc
rrcFilter = rcosdesign(rolloff,span,sps,'normal');
txPulse = upfirdn(modData, rrcFilter, sps)*sqrt(sps);% power = 1
% eyediagram(txPulse(sps*3+1:2000), sps)
sig = plotEye(txPulse(sps*3:2001), sps, sps/2+2);
hold
plot([0 0], [-1.7 1.7],'Color',[0.8500    0.3250    0.0980])
text(0-0.02, -1.8, '(1)' )
plot([0.125 0.125], [-1.7 1.7],'Color',[0.9290    0.6940    0.1250])
text(0.125-0.02, -1.8, '(2)' )
plot([0.25 0.25], [-1.7 1.7],'Color',[0.4940    0.1840    0.5560])
text(0.25-0.02, -1.8, '(3)' )
plot([0.375 0.375], [-1.7 1.7],'Color',[0.4660    0.6740    0.1880])
text(0.375-0.02, -1.8, '(4)' )
myplotsetting
hold off
% we will need histogram
nbins = 18;
ylimit = 0.6;
figure
subplot(2,2,1)
h=histogram(sig(17,:),nbins,'Normalization','probability');
h.FaceColor = [0.8500    0.3250    0.0980];
ylim([0 ylimit])
ylabel('Probability')
xlabel('(1)')
% myplotsetting
subplot(2,2,2)
h=histogram(sig(21,:),nbins,'Normalization','probability');
h.FaceColor = [0.9290    0.6940    0.1250];
ylim([0 ylimit])
ylabel('Probability')
% myplotsetting
xlabel('(2)')
subplot(2,2,3)
myplotsetting
h=histogram(sig(25,:),nbins,'Normalization','probability');
h.FaceColor = [0.4940    0.1840    0.5560];
ylim([0 ylimit])
ylabel('Probability')
xlabel('(3)')
% myplotsetting
subplot(2,2,4)
h=histogram(sig(29,:),nbins,'Normalization','probability');
h.FaceColor = [0.4660    0.6740    0.1880];
ylim([0 ylimit])
ylabel('Probability')
xlabel('(4)')
% myplotsetting
% figure
% histogram(sig(17,:),nbins,'Normalization','probability');
% ylim([0 ylimit])
% ylabel('Probability')
% hold on
% histogram(sig(21,:),nbins,'Normalization','probability');
% histogram(sig(25,:),nbins,'Normalization','probability');
% histogram(sig(29,:),nbins,'Normalization','probability');
% hold off
