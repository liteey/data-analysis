
% trim the vectors so they are the same length
trimLen = min([length(cam1x), length(cam2x), length(cam3x)]);
trimCam1x = cam1x(1:trimLen);
trimCam1y = cam1y(1:trimLen);
trimCam2x = cam2x(1:trimLen);
trimCam2y = cam2y(1:trimLen);
trimCam3x = cam3x(1:trimLen);
trimCam3y = cam3y(1:trimLen);

% Perform the Principal Component Analysis
firstMat = [trimCam1x trimCam1y trimCam2x trimCam2y trimCam3x trimCam3y]';
% SVD of X will give you one singular value if it is perfect data

[U,S,V] = svd(firstMat, 'econ'); 
sig = diag(S);

figure(1)
% subplot(1,2,1)
plot(sig.^2/sum(sig.^2),'ko','Linewidth',1)
xlabel('Sigma Index'); ylabel('Energy')
title('Sigma Energies of Ideal Case')

% prin_comp = U' * firstMat;
% plot(1:trimLen, prin_comp(1:3,:),'Linewidth',1)
% xlabel('Time (frames)'); ylabel('Displacement (pixels)')
% legend('1st Principal Component', '2nd Principal Component', ...
%        '3rd Principal Component', 'location', 'northeast')
% title('Case 1: First Three Principal Component ProjectionS')
%% Test 2: Noisy Case
clear; clc;
load('cam1_2.mat')
load('cam2_2.mat')
load('cam3_2.mat')

% implay(vidFrames1_2)
cropVidFrames1_2 = vidFrames1_2(:,300:450,:,:);
% implay(cropVidFrames1_2)
cropVidFrames2_2 = vidFrames2_2(:,180:450,:,:);
% implay(cropVidFrames2_2)
cropVidFrames3_2 = vidFrames3_2(175:350,:,:,:);
% implay(cropVidFrames3_2)

numFrames_1 = size(cropVidFrames1_2,4);
cam1 = zeros(numFrames_1,2);
for i = 1:numFrames_1
    t1 = cropVidFrames1_2(:,:,1,i);
    [Max,Ind] = max(t1(:));
    [x1,y1] = ind2sub([size(t1,1), size(t1,2)], Ind);
    cam1(i,:) = [x1,y1];
end

numFrames_2 = size(cropVidFrames2_2,4);
cam2 = zeros(numFrames_2,2);
for i = 1:numFrames_2
    t2 = cropVidFrames2_2(:,:,1,i);
    [Max,Ind] = max(t2(:));
    [x1,y1] = ind2sub([size(t2,1), size(t2,2)], Ind);
    cam2(i,:) = [x1,y1];
end

