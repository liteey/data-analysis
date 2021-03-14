clear; close all; clc;

[images, labels] = mnist_parse('train-images.idx3-ubyte', 'train-labels.idx1-ubyte');
col_images = zeros(784,60000);
for k=1:60000
%    subplot(3,3,k)
   num1 = reshape(images(:,:,k),784,1);
   col_images(:,k) = num1;
end

[m,n] = size(col_images);
for k=1:n
    col_images(:,k) = col_images(:,k) - mean(col_images(:,k));
% mn = mean(col_images,1);
% col_images = col_images - mn;

end

[U,S,V] = svd(col_images, 'econ');


%% plot singular values
sig = diag(S);
% 
% figure(1)
% plot(sig.^2/sum(sig.^2),'ko','Linewidth',1)
% xlabel('Sigma Index'); ylabel('Energy')
% title('Sigma Energies')

energy = 0;
iters = 1;
while energy < 0.9
    energy = energy + (sig(iters).^2/sum(sig.^2));
    iters = iters + 1;
    
end
iters = iters - 1;

% seems like we need 64 modes to grab 90% of the energy
%% project onto PCA modes
% sigma times v
% first col of v, one is pos one is neg, check different cols of v
% inform which is hard to separate, which is easy to separate

[unique_groups, ~, group_idx] = unique(labels);
num_groups= size(unique_groups,1);

fullproj = U'*col_images;
%% Sorting based on label
ind0 = 1;
zeroInd = zeros();
ind1 = 1;
oneInd = zeros();
ind2 = 1;
twoInd = zeros();
ind3 = 1;
threeInd = zeros();
ind4 = 1;
fourInd = zeros();
ind5 = 1;
fiveInd = zeros();
ind6 = 1;
sixInd = zeros();
ind7 = 1;
sevenInd = zeros();
ind8 = 1;
eightInd = zeros();
ind9 = 1;
nineInd = zeros();

for k=1:n
    if(unique_groups(group_idx(k)) == 0)
        zeroInd(ind0) = k;
        ind0 = ind0 + 1;
    elseif(unique_groups(group_idx(k)) == 1)
        oneInd(ind1) = k;
        ind1 = ind1 + 1;
    elseif(unique_groups(group_idx(k)) == 2)
        twoInd(ind2) = k;
        ind2 = ind2 + 1;
    elseif(unique_groups(group_idx(k)) == 3)
        threeInd(ind3) = k;
        ind3 = ind3 + 1;
    elseif(unique_groups(group_idx(k)) == 4)
        fourInd(ind4) = k;
        ind4 = ind4 + 1;
    elseif(unique_groups(group_idx(k)) == 5)
        fiveInd(ind5) = k;
        ind5 = ind5 + 1;
    elseif(unique_groups(group_idx(k)) == 6)
        sixInd(ind6) = k;
        ind6 = ind6 + 1;
    elseif(unique_groups(group_idx(k)) == 7)
        sevenInd(ind7) = k;
        ind7 = ind7 + 1;
    elseif(unique_groups(group_idx(k)) == 8)
        eightInd(ind8) = k;
        ind8 = ind8 + 1;
    elseif(unique_groups(group_idx(k)) == 9)
        nineInd(ind9) = k;
        ind9 = ind9 + 1;
    end
end

%%
zs = fullproj(1:iters,zeroInd);
ones = fullproj(1:iters,oneInd);
twos = fullproj(1:iters,twoInd);
threes = fullproj(1:iters,threeInd);
fours = fullproj(1:iters,fourInd);
fives = fullproj(1:iters,fiveInd);
sixs = fullproj(1:iters,sixInd);
sevens = fullproj(1:iters,sevenInd);
eights = fullproj(1:iters,eightInd);
nines = fullproj(1:iters,nineInd);

%%
hold on
plot3(  ones(3,:),   ones(2,:),   ones(5,:),'ko')
plot3(    zs(3,:),     zs(2,:),     zs(5,:),'ro')
plot3(  twos(3,:),   twos(2,:),   twos(5,:),'go')
plot3(threes(3,:), threes(2,:), threes(5,:),'yo')
plot3( fours(3,:),  fours(2,:),  fours(5,:),'co')
plot3( fives(3,:),  fives(2,:),  fives(5,:),'bo')
plot3(  sixs(3,:),   sixs(2,:),   sixs(5,:),'o', 'Color', [.61 .51 .74])
plot3(sevens(3,:), sevens(2,:), sevens(5,:),'wo')
plot3(eights(3,:), eights(2,:), eights(5,:),'o', 'Color', [.31 .41 .92])
plot3( nines(3,:),  nines(2,:),  nines(5,:),'o', 'Color', [.87 .65 .11])
% proj = S*V(:,1)';


% seven and nine stick together
% pro
%% Build LDA for two digit classification
numSev = size(sevenInd,2);
numNin = size(nineInd,2);
meanSev = mean(sevens,2);
meanNin = mean(nines,2);

wn = 0; % within variances
for k = 1:numSev
   wn = wn + (sevens(:,k)-meanSev)*(sevens(:,k)-meanSev)'; 
end

for k = 1:numNin
    wn = wn + (nines(:,k) - meanNin)*(nines(:,k)-meanNin)';
end

bw = (meanSev-meanNin)*(meanSev-meanNin)';

% finding proj using LDA
[V2, D] = eig(bw, wn); 
[lambda, ind] = max(abs(diag(D)));
w = V2(:,ind);
w = w/norm(w,2);

projSev = w'*sevens;
projNin = w'*nines;

if mean(projSev) > mean(projNin)
    w = -w;
    projSev = -projSev;
    projNin = -projNin;
end

% nines below 0
% sevens above 0
figure
plot(projSev,sevens(length(projSev)),'ko')
hold on
plot(projNin,nines(length(projNin)),'bo')

% there is pretty big overlap between the values
%%
% Find threshold
sortSev = sort(projSev);
sortNin = sort(projNin);

t1 = length(sortSev);
t2 = 1;
while sortSev(t1) > sortNin(t2)
    t1 = t1 - 1;
    t2 = t2 + 1;
end
threshold = (sortSev(t1) + sortNin(t2))/2;

%% histogram of seven and nine values

figure
subplot(1,2,1)
histogram(sortSev,30); hold on, plot([threshold threshold], [0 1200], 'r')
title('Seven')
subplot(1,2,2)
histogram(sortNin,30); hold on, plot([threshold threshold], [0 1200], 'r')
title('Nine')
%% Test the data
[images1, labels1] = mnist_parse('t10k-images.idx3-ubyte', 't10k-labels.idx1-ubyte');

% TestSet = zeros();
TestSet = zeros(784, 10000);
for k=1:10000
   num1 = reshape(images1(:,:,k),784,1);
   TestSet(:,k) = num1;
end

[g, ~, gidx] = unique(labels1);
sevenIndTest = zeros();
nineIndTest = zeros();
ind7Test = 1;
ind9Test = 1;

for k=1:10000
    if(g(gidx(k)) == 7)
        sevenIndTest(ind7Test) = k;
        ind7Test = ind7Test + 1;
    elseif(g(gidx(k)) == 9)
        nineIndTest(ind9Test) = k;
        ind9Test = ind9Test + 1;
    end
end

TestMat = U'*TestSet; % PCA projection
TestMat = TestMat(1:iters,:); % go up to 90% energy

sevensTest = TestMat(1:iters,sevenIndTest);
ninesTest = TestMat(1:iters,nineIndTest);

testcase = [sevensTest, ninesTest];

LabelsVec = zeros(1,size(testcase,2));
for k = size(sevensTest,2):size(testcase,2)
   LabelsVec(k) = 1; 
end

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
while sortSevNin(t1) > sortFou(t2)
    t1 = t1 - 1;
    t2 = t2 + 1;
end
threshold = (sortSevNin(t1) + sortFou(t2))/2;

figure
subplot(1,2,1)
histogram(sortSevNin,30); hold on, plot([threshold threshold], [0 2000], 'r')
title('Seven Nine')
subplot(1,2,2)
histogram(sortFou,30); hold on, plot([threshold threshold], [0 1200], 'r')
title('Four')

% sv - maximize this quantity - max the mean between the two datasets -
% dictates how much overlap there is
% want no relation between class 1 and class 2 - correlation in each class
% sw - variance matrix - variance in each class
%% Test the data
[images1, labels1] = mnist_parse('t10k-images.idx3-ubyte', 't10k-labels.idx1-ubyte');

% TestSet = zeros();
TestSet = zeros(784, 10000);
for k=1:10000
   num1 = reshape(images1(:,:,k),784,1);
   TestSet(:,k) = num1;
end

[g, ~, gidx] = unique(labels1);
sevenIndTest = zeros();
nineIndTest = zeros();
fourIndTest = zeros();
ind7Test = 1;
ind9Test = 1;
ind4Test = 1;

for k=1:10000
    if(g(gidx(k)) == 7)
        sevenIndTest(ind7Test) = k;
        ind7Test = ind7Test + 1;
    elseif(g(gidx(k)) == 9)
        nineIndTest(ind9Test) = k;
        ind9Test = ind9Test + 1;
    elseif(g(gidx(k)) == 4)
        fourIndTest(ind4Test) = k;
        ind4Test = ind4Test + 1;
    end
end


TestMat = U'*TestSet; % PCA projection
TestMat = TestMat(1:iters,:); % go up to 90% energy

sevensTest = TestMat(1:iters,sevenIndTest);
ninesTest = TestMat(1:iters,nineIndTest);
foursTest = TestMat(1:iters,fourIndTest);

testcase = [sevensTest, ninesTest, foursTest];

% seven and nine is 0, four is 1
LabelsVec = zeros(1,size(testcase,2));
for k = size(sevensTest,2)+size(ninesTest,2):size(testcase,2)
   LabelsVec(k) = 1; 
end

TestNum = size(testcase,2);
pval = w'*testcase;

ResVec = (pval > threshold);
err = abs(ResVec - LabelsVec);
errNum = sum(err);
sucRate = 1 - errNum/TestNum; % 0.9506




%% Classification Tree
clear; close all; clc;
% classification tree on fisheriris data
% load fisheriris;
% tree=fitctree(meas,species,'MaxNumSplits',3,'CrossVal','on');
% view(tree.Trained{1},'Mode','graph');
% classError = kfoldLoss(tree)

[images, labels] = mnist_parse('train-images.idx3-ubyte', 'train-labels.idx1-ubyte');
[test, labels1] = mnist_parse('t10k-images.idx3-ubyte', 't10k-labels.idx1-ubyte');

col_images = zeros(784,60000);
for k=1:60000
%    subplot(3,3,k)
   num1 = reshape(images(:,:,k),784,1);
   col_images(:,k) = num1;
end

small = col_images(:,1:60000);
smalllab = labels(1:60000,1);

tree = fitctree(small',smalllab,'MaxNumSplits',10,'CrossVal','on');
view(tree.Trained{1},'Mode','graph');
classError = kfoldLoss(tree);




%% Using SVM
[images, labels] = mnist_parse('train-images.idx3-ubyte', 'train-labels.idx1-ubyte');
[test, labels1] = mnist_parse('t10k-images.idx3-ubyte', 't10k-labels.idx1-ubyte');
col_images = zeros(784,60000);
for k=1:60000
%    subplot(3,3,k)
   num1 = reshape(images(:,:,k),784,1);
   col_images(:,k) = num1;
end

small = col_images(:,1:1000);
smalllab = labels(1:1000,1);


test_images = zeros(784,10000);
for k=1:10000
%    subplot(3,3,k)
   num1 = reshape(test(:,:,k),784,1);
   test_images(:,k) = num1;
end

smalltest = test_images(:,1:1000);

% SVM classifier with training data, labels and test set
Mdl = fitcecoc(small',smalllab);
test_labels = predict(Mdl,smalltest');

err_num = 0;
for k = 1:size(test_labels,1)
   if (test_labels(k) - labels(k) ~= 0)
      err_num = err_num + 1;
   end
end

total_err = err_num;
sucrate = 1 - err_num/size(test_labels,1);

%%
[images, labels] = mnist_parse('train-images.idx3-ubyte', 'train-labels.idx1-ubyte');



%%
[images1, labels1] = mnist_parse('t10k-images.idx3-ubyte', 't10k-labels.idx1-ubyte');
