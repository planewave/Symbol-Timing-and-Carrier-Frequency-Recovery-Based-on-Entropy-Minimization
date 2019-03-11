function entp = simpEntpBat(sig,r, gpu)
narginchk(2,3);
if nargin == 2
    gpu = false;
end


len = length(sig);
sig = sig(:);
sig = [real(sig) imag(sig)];

if gpu == 1
    sig = gpuArray(single(sig));
end

dis = pdist(sig, 'euclidean');
entp = sum(dis>r);
entp = entp/(len*(len-1)/2);