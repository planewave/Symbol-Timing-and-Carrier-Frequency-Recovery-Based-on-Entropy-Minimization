function sig = plotEye(sig,sps,offset)
%----------------------------------------------
%Function to plot eye diagram of inputSignal
% real component only
if(nargin < 3),offset=1;end
sig = real(sig(offset:end));
nRep = floor(length(sig)/sps); % 必须是整数!
sig(sps*nRep+1) = NaN;
sig = sig(1:sps*nRep+1); % 这里需要多一个数据点

eyeX = (-sps/2:sps/2)'/sps;
% eyeX = (-sps/2:sps/2)';

sig_rs = reshape(sig(1:end-1),[sps,nRep]);
sig_fix = [sig_rs(1,:),sig(end)];
sig = [sig_rs;sig_fix(2:end)];
% plot(eyeX,sig,'-','Color',[0 0.4470 0.7410]);
plot(eyeX,sig,'-b'); % for black

xlabel('Normalized symbol period')
ylabel('Amplitude')
% out = 1;