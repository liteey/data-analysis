
TestNum = size(testcase,2);
pval = w'*testcase;

ResVec = (pval > threshold);
err = abs(ResVec - LabelsVec);
errNum = sum(err);
sucRate = 1 - errNum/TestNum;

% still have quite a bit that is mis classified
%% Build LDA for Classify three digits
SevsNins = [sevens,nines];
SevNinInd = [sevenInd, nineInd];
numSevNin = numSev + numNin;

numFou = size(fourInd,2);
meanFou = mean(fours,2);
meanSevNin = mean(SevsNins,2);

wn = 0; % within variances
for k = 1:numSevNin
   wn = wn + (SevsNins(:,k)-meanSevNin)*(SevsNins(:,k)-meanSevNin)'; 
end

for k = 1:numFou
    wn = wn + (fours(:,k) - meanFou)*(fours(:,k)-meanFou)';
end

bw = (meanSevNin-meanFou)*(meanSevNin-meanFou)';

% LDA
[V2, D] = eig(bw,wn);
[lambda, ind] = max(abs(diag(D)));
w = V2(:,ind);
w = w/norm(w,2);

projSevNin = w'*SevsNins;
projFou = w'*fours;

if mean(projSevNin) > mean(projFou)
    w = -w;
    projSevNin = -projSevNin;
    projFou = -projFou;
end

figure
plot(projSevNin,SevsNins(length(projSevNin)),'ko')
hold on
plot(projFou,fours(length(projFou)),'bo')

% a lot more overlap

%% 

sortSevNin = sort(projSevNin);
sortFou = sort(projFou);

t1 = length(sortSevNin);
t2 = 1;

