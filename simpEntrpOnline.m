function offset = simpEntrpOnline(sig,r)

window = 100;
colNum = 5;
persistent sigMat distMat entArray;

if isempty(sigMat)
    sigMat = zeros(window,colNum);
    distMat = zeros(window,colNum);
    entArray = zeros(1,colNum);
end

sig = sig(:).';
sigMat = [sig;sigMat(1:end-1,:)];
distMat = [dist(sigMat,r);distMat(1:end-1,:)];
entArray = entArray-distMat(end,:)+distMat(1,:);

if distMat(end,2) == 0 || distMat(end,4) == 0
    % the timing offset will not change
    offset = 0;
elseif entArray(3)>entArray(2)
    offset = -1;
    sigMat = right(sigMat);
    distMat = right(distMat);
    entArray = right(entArray);
elseif entArray(3)>entArray(4)
    offset = 1;
    sigMat = left(sigMat);
    distMat = left(distMat);
    entArray = left(entArray);
else
    offset = 0;
end
end
%% assitant functions
function d = dist(sigArray,r) % Euclid distance
dif = abs(bsxfun(@minus,sigArray(2:end,:),sigArray(1,:))); 
d = sum(dif>r);
end

function newMat = right(mat)
    newMat = zeros(size(mat));
	newMat(:,2:5) = mat(:,1:4);
end

function newMat = left(mat)
    newMat = zeros(size(mat));
	newMat(:,1:4) = mat(:,2:5);
end