function entp = bmre(sig, r, r2, gpu)
% BMRE estimation in batch
% each data set must be in one column 
% r entropy threshold
narginchk(3, 4);
if nargin == 3
    gpu = false;
end

if gpu == true
    sig = gpuArray(single(sig));
end

[row, col] = size(sig);
sig2 = real(zeros(row,col*2,'like',sig));
% 'like' to make it same cpu or gpu array as sig
% real() and the following are to make it possible
% to run pdist()
sig2(:,1:2:end) = real(sig);
sig2(:,2:2:end) = imag(sig);

sigCell = mat2cell(sig2,row,2*ones(1,col));
% clear sig sig2
% sig = sig(abs(sig)>r2);
% sig = sig(:); % convert to column
% dis = pdist(sig, 'euclidean'); 
amp = @(x) x(x>r2);
sigCell = cellfun(@(x) amp(x), sigCell,'UniformOutput', false);

dist = cellfun(@(x) pdist(x), sigCell, ...
    'UniformOutput', false);
% entp = gether(dist);

entp = cellfun(@(x) 1-sum(x<r)/(row*(row-1)/2),dist);