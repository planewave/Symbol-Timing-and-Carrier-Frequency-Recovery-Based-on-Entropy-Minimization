function out = plotfft(sig,fs,order,nfft)

narginchk(2,4);
len = length(sig);
if nargin == 2
    order = 1;
    nfft = len;
elseif nargin == 3
    nfft = len;
    
    
    if rem(len,2) == 1
        len2 = (len-1)/2;
    else
        len2 = (len-2)/2;
        sig = sig(1:end-1);
        nfft = nfft-1;
    end
    figure
    plot((-len2:len2)/len*fs,(20*log10(abs(fftshift((fft(sig.^order,nfft)))))));
    myplotsetting
    
    elseif nargin == 4
        if rem(nfft,2) == 1
            len2 = (nfft-1)/2;
        else
            len2 = nfft/2;
%             sig = sig(1:end-1);
%             nfft = nfft-1;
        end
        figure
        plot((-len2:len2)/nfft*fs,(20*log10(abs(fftshift((fft(sig.^order,nfft)))))));
        myplotsetting
end
out = 1;