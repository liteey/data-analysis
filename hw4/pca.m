for k=1:4
    subplot(2,2,k)
    ut1 = reshape(U(:,k),32,32);
    ut2 = rescale(ut1);
    imshow(ut2)
end