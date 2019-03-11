function entp = simpEntpEng(sig,r,r2,gpu)
% entropy estimation with energy threshold
% r entropy threshold, r2 energy threshold
narginchk(3,4);
    if nargin == 3
        gpu = false;
    end
len = length(sig);
sig = sig(abs(sig)>r2);
sig = sig(:); % convert to column
sig = [real(sig) imag(sig)];
if gpu == true
    sig = gpuArray(single(sig));
end
dis = pdist(sig, 'euclidean'); 
entp = 1 - sum(dis<r)/(len*(len-1)/2);