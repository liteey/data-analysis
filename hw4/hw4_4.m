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
