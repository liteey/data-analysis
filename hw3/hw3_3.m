numFrames_3 = size(cropVidFrames3_2,4);
cam3 = zeros(numFrames_3,2);
for i = 1:numFrames_3
    t3 = cropVidFrames3_2(:,:,1,i);
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

secondMat = [trimCam1x trimCam1y trimCam2x trimCam2y trimCam3x trimCam3y]';
[U,S,V] = svd(secondMat, 'econ'); 
sig = diag(S);
figure(2)
plot(sig.^2/sum(sig.^2),'ko','Linewidth',1)
xlabel('Sigma Index'); ylabel('Energy')
title('Sigma Energies of Noisy Case')

% prin_comp = U' * secondMat;
% plot(1:trimLen, prin_comp(1:3,:),'Linewidth',1)
% xlabel('Time (frames)'); ylabel('Displacement (pixels)')
% legend('1st Principal Component', '2nd Principal Component', ...
%        '3rd Principal Component', 'location', 'northeast')
% title('Case 2: First Three Principal Component ProjectionS')

%% Test 3: Horizontal Displacement
clear; clc;
load('cam1_3.mat')
load('cam2_3.mat')
load('cam3_3.mat')


cropVidFrames1_3 = vidFrames1_3(:,250:400,:,:);
% implay(cropVidFrames1_3)
cropVidFrames2_3 = vidFrames2_3(:,180:450,:,:);
% implay(cropVidFrames2_3)
cropVidFrames3_3 = vidFrames3_3(150:350,:,:,:);
% implay(cropVidFrames3_3)

numFrames_1 = size(cropVidFrames1_3,4);
cam1 = zeros(numFrames_1,2);
for i = 1:numFrames_1
    t1 = cropVidFrames1_3(:,:,1,i);
    [Max,Ind] = max(t1(:));
    [x1,y1] = ind2sub([size(t1,1), size(t1,2)], Ind);
    cam1(i,:) = [x1,y1];
end
