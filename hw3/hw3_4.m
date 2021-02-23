
numFrames_2 = size(cropVidFrames2_3,4);
cam2 = zeros(numFrames_2,2);
for i = 1:numFrames_2
    t2 = cropVidFrames2_3(:,:,1,i);
    [Max,Ind] = max(t2(:));
    [x1,y1] = ind2sub([size(t2,1), size(t2,2)], Ind);
    cam2(i,:) = [x1,y1];
end

numFrames_3 = size(cropVidFrames3_3,4);
cam3 = zeros(numFrames_3,2);
for i = 1:numFrames_3
    t3 = cropVidFrames3_3(:,:,1,i);
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

thirdMat = [trimCam1x trimCam1y trimCam2x trimCam2y trimCam3x trimCam3y]';
[U,S,V] = svd(thirdMat, 'econ'); 
sig = diag(S);
figure(3)
plot(sig.^2/sum(sig.^2),'ko','Linewidth',1)
title('Sigma Energies of Horizontal Movement Case')
xlabel('Sigma Index'); ylabel('Energy')

% prin_comp = U' * thirdMat;
% plot(1:trimLen, prin_comp(3:5,:),'Linewidth',1)
% xlabel('Time (frames)'); ylabel('Displacement (pixels)')
% % legend('1st Principal Component', '2nd Principal Component', ...
% %        '3rd Principal Component', '4th Principal Component', ...
% %        '5th Principal Component', 'location', 'southeast')
%  legend('3rd Principal Component', '4th Principal Component', ...
%         '5th Principal Component', 'location', 'southeast')
% title('Case 3: Three to Five Principal Component ProjectionS')

%% Part 4: Horizontal Displacement and Rotation
%clear; close all; clc;
clear; clc;
load('cam1_4.mat')
load('cam2_4.mat')
load('cam3_4.mat')

cropVidFrames1_4 = vidFrames1_4(:,300:500,:,:);
% implay(cropVidFrames1_4)
