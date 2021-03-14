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


