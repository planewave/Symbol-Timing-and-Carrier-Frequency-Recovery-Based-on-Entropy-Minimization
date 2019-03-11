function [y, err, phi, lambda] = carrSync(x, ref)
    % carrier sync for BPSK, QPSK
    narginchk(1, 2);
    len = length(x);
    y = x;
    err = zeros(1,len);
    Bn = 0.08; % NormalizedLoopBandwidth
    xi = 0.707; % DampingFactor
    Kp = 2; % for BPSK QPSK
    K0 = 1; % SamplesPerSymbol
    phi = zeros(1,len);
    lambda = zeros(1,len);
    theta = Bn/(xi+1/(4*xi));
    d = 1+2*xi*theta+theta^2;
    g1 = 4*(theta^2/d)/(Kp*K0);
    gp = 4*xi*(theta/d)/(Kp*K0);
    for n  = 2:len
        % phase shift
        y(n) = x(n)*exp(-1i*lambda(n-1));
        % PLL there has to be both positive and negtive 
%         refDec = exp(1j*(pskdemod(y(n),4,pi/4)*pi/2+pi/4)); 
%         err(n) = angle(y(n)*refDec'); % DDA
        if nargin == 2
            err(n) = angle(y(n)*ref(n)'); % DA
        else
            err(n) = sign(real(y(n)))*imag(y(n)) ...
            -sign(imag(y(n)))*real(y(n)); % NDA
        end
        

%         
        
%         err(n) = sign(real(x(n)))*imag(x(n)); % uncomment for BPSK
        phi(n) = g1*err(n) + phi(n-1);
        % DDS
        lambda(n) = (gp*err(n-1)+phi(n-1))+lambda(n-1);
    end    

        
end