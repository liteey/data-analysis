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