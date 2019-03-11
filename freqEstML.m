function [f_est, Data_comp, pha_out]= freqEstML(offsetData,fs,order,pha_in)
% ML based CFO global search
% +/- 0.5 Hz search range. Note that the raised signal reduce the range
% for order = 4, actual range is +/- 0.125 Hz
% fs sampling rate; pha_in initial phase
    narginchk(3,4);
    if nargin == 3
        pha_in = 0;
    end
    raisedSignal = offsetData(:).^ order;
    len = length(raisedSignal);
    est_len = 101;
    est_start  = -0.51; 

    eng = zeros(1,est_len);
    for n = 1:est_len
        f_tr = n/100 + est_start;
        x_tr = raisedSignal.*exp(-1j*2*pi*f_tr*(0:len-1)/fs).';
%         entrp(n) = simpEntrpBat(x_tr,r);
        eng(n) = abs(sum(x_tr));
    end
%     plot(-0.5:0.01:0.5,eng)
    [~,ind] = max(eng);
    f_est = est_start+ind/100;

    for n = 1:21  
        f_tr = n/1000+f_est-0.011;
        x_tr = raisedSignal.*exp(-1j*2*pi*f_tr*(0:len-1)/fs).';
        eng(n) = abs(sum(x_tr));
    end
    [~,ind] = max(eng);
    f_est = ind/1000+f_est-0.011;
    
    f_est = f_est/order;
    Data_comp = offsetData.*exp(-1j*(2*pi*(f_est)*(0:len-1)/fs+pha_in)).';
    [Data_comp, pha] = phase_restore(Data_comp);
    pha_out = 2*pi*(f_est)*(len)/fs+pha_in+pha;
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