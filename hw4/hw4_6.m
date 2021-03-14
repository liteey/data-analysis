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

