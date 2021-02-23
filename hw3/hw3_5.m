cropVidFrames2_4 = vidFrames2_4(:,180:450,:,:);
% implay(cropVidFrames2_4)
cropVidFrames3_4 = vidFrames3_4(150:350,:,:,:);
% implay(cropVidFrames3_4)

numFrames_1 = size(cropVidFrames1_4,4);
cam1 = zeros(numFrames_1,2);
for i = 1:numFrames_1
    t1 = cropVidFrames1_4(:,:,1,i);
    [Max,Ind] = max(t1(:));
    [x1,y1] = ind2sub([size(t1,1), size(t1,2)], Ind);
    cam1(i,:) = [x1,y1];
end

numFrames_2 = size(cropVidFrames2_4,4);
cam2 = zeros(numFrames_2,2);
for i = 1:numFrames_2
    t2 = cropVidFrames2_4(:,:,1,i);
    [Max,Ind] = max(t2(:));
    [x1,y1] = ind2sub([size(t2,1), size(t2,2)], Ind);
    cam2(i,:) = [x1,y1];
end

numFrames_3 = size(cropVidFrames3_4,4);
cam3 = zeros(numFrames_3,2);
for i = 1:numFrames_3
    t3 = cropVidFrames3_4(:,:,1,i);
    [Max,Ind] = max(t3(:));
    [x1,y1] = ind2sub([size(t3,1), size(t3,2)], Ind);
    cam3(i,:) = [x1,y1];
end

cam1x = cam1(:,1) - mean(cam1(:,1));
cam1y = cam1(:,2) - mean(cam1(:,2));
cam2x = cam2(:,1) - mean(cam2(:,1));
cam2y = cam2(:,2) - mean(cam2(:,2));
cam3x = cam3(:,1) - mean(cam3(:,1));
cam3y = cam3(:,2) - mean(cam3(:,2));

trimLen = min([length(cam1x), length(cam2x), length(cam3x)]);
trimCam1x = cam1x(1:trimLen);
trimCam1y = cam1y(1:trimLen);
trimCam2x = cam2x(1:trimLen);
trimCam2y = cam2y(1:trimLen);
trimCam3x = cam3x(1:trimLen);
trimCam3y = cam3y(1:trimLen);


fourthMat = [trimCam1x trimCam1y trimCam2x trimCam2y trimCam3x trimCam3y]';
[U,S,V] = svd(fourthMat, 'econ'); 
sig = diag(S);
figure(4)
% plot(sig.^2/sum(sig.^2),'ko','Linewidth',1)
% title('Sigma Energies of Horizontal Movement and Rotation Case')
% xlabel('Sigma Index'); ylabel('Energy')

