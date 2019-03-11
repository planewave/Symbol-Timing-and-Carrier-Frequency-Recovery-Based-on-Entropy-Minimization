function [f_est, Data_comp, pha_out]= freqEstEntp(offsetData, fs, r, pha_in, gpu)
% EM based CFO search. range: 0.5 Hz, accuracy: 0.001 Hz
% fs: sampling rate. pha_in: initial phase
    narginchk(3, 5);
    if nargin == 3
        pha_in = 0;
        gpu = 0;
    elseif nargin == 4
        gpu = 0;
    end
    offsetData = offsetData(:);
    len = length(offsetData);
    est_len = 101;
    est_start  = -0.51;

%     entrp = zeros(1,est_len);
    parfor n = 1:est_len
        f_tr = n/100 + est_start;
        x_tr = offsetData.*exp(-1j*2*pi*f_tr*(0:len-1)/fs).';
        entp(n) = simpEntpBat(x_tr, r, gpu);
    end
%     plot(-0.5:0.01:0.5,entp)
    [~,ind] = min(entp);
    f_est = est_start+ind/100;
    
%     entp = [];
    parfor n = 1:21 
        f_tr = n/1000+f_est-0.011;
        x_tr = offsetData.*exp(-1j*2*pi*f_tr*(0:len-1)/fs).';
        entp2(n) = simpEntpBat(x_tr, r, gpu);
    end
    [~,ind] = min(entp2);
    f_est = ind/1000+f_est-0.011;

    Data_comp = offsetData.*exp(-1j*(2*pi*(f_est)*(0:len-1)/fs+pha_in)).';
    [Data_comp, pha] = phase_restore(Data_comp);
    pha_out = 2*pi*(f_est)*(len)/fs+pha_in+pha;
end

function entp = simpEntpBat(sig, r, gpu)
%     old algorithm commented 
%     len = length(sig);
%     entrp = 0;
% 
%     for n = 1:len-1
%         dif = abs(sig(n)-sig(n+1:end));% Euclidean distance
%         entrp = entrp+sum(dif>r);
%     end
%     entrp = entrp/(len*(len-1)/2);
    len = length(sig);
    sig = sig(:);
    sig = [real(sig) imag(sig)];

    if gpu == 1
        sig = gpuArray(single(sig));
    end

    dis = pdist(sig, 'euclidean');
    entp = sum(dis>r);
    entp = entp/(len*(len-1)/2);
end

function [sig, pha] = phase_restore(x)
% phase recovery
% Rice, Michael. Digital Communications: A Discrete-Time Approach. Upper Saddle River, NJ: Prentice Hall, 2009, pp. 359?393.
    sig = x;
    pha = 0;
    err = 1;
    while abs(err) > 0.05
        err = mean((sign(real(sig)).*imag(sig)-(sign(imag(sig)).*real(sig))));
        sig = sig*exp(-1j*err*0.1);
        pha = pha+err*0.1;
    end
end