%% Test 1: Ideal Case
clear; close all; clc;

load('cam1_1.mat')
% implay(vidFrames1_1)
load('cam2_1.mat')
% implay(vidFrames2_1)
load('cam3_1.mat')
% implay(vidFrames3_1)


cropVidFrames1_1 = vidFrames1_1(:,275:400,:,:);
% implay(cropVidFrames1_1)
cropVidFrames2_1 = vidFrames2_1(:,200:370,:,:);
% implay(cropVidFrames2_1)
cropVidFrames3_1 = vidFrames3_1(210:350,:,:,:);
% implay(cropVidFrames3_1)



numFrames_1 = size(cropVidFrames1_1,4);
cam1 = zeros(numFrames_1,2);
for i = 1:numFrames_1
    t1 = cropVidFrames1_1(:,:,1,i);
    [Max,Ind] = max(t1(:));
    [x1,y1] = ind2sub([size(t1,1), size(t1,2)], Ind);
    cam1(i,:) = [x1,y1];
end

numFrames_2 = size(cropVidFrames2_1,4);
cam2 = zeros(numFrames_2,2);
for i = 1:numFrames_2
    t2 = cropVidFrames2_1(:,:,1,i);
    [Max,Ind] = max(t2(:));
    [x1,y1] = ind2sub([size(t2,1), size(t2,2)], Ind);
    cam2(i,:) = [x1,y1];
end

numFrames_3 = size(cropVidFrames3_1,4);
cam3 = zeros(numFrames_3,2);
for i = 1:numFrames_3
    t3 = cropVidFrames3_1(:,:,1,i);
    [Max,Ind] = max(t3(:));
    [x1,y1] = ind2sub([size(t3,1), size(t3,2)], Ind);
    cam3(i,:) = [x1,y1];
end

% Set mean of the results equal to 0
cam1x = cam1(:,1) - mean(cam1(:,1));
cam1y = cam1(:,2) - mean(cam1(:,2));
cam2x = cam2(:,1) - mean(cam2(:,1));
cam2y = cam2(:,2) - mean(cam2(:,2));
cam3x = cam3(:,1) - mean(cam3(:,1));
cam3y = cam3(:,2) - mean(cam3(:,2));

% figure(1)
% plot(cam1x)
% hold on
% plot(cam1y)
% 
% figure(2)
% plot(cam2x)
% hold on
% plot(cam2y)
% 
% figure(3)
% plot(cam3x)
% hold on
% plot(cam3y)

% trim the vectors so they are the same length
trimLen = min([length(cam1x), length(cam2x), length(cam3x)]);
trimCam1x = cam1x(1:trimLen);
trimCam1y = cam1y(1:trimLen);
trimCam2x = cam2x(1:trimLen);
trimCam2y = cam2y(1:trimLen);
trimCam3x = cam3x(1:trimLen);
trimCam3y = cam3y(1:trimLen);

% figure(1)
% plot(trimCam1x)
% hold on
% plot(trimCam1y)
% 
% figure(2)
% plot(trimCam2x)
% hold on
% plot(trimCam2y)
% 
% figure(3)
% plot(trimCam3x)
% hold on
% plot(trimCam3y)


% Perform the Principal Component Analysis
firstMat = [trimCam1x trimCam1y trimCam2x trimCam2y trimCam3x trimCam3y]';
% SVD of X will give you one singular value if it is perfect data

[U,S,V] = svd(firstMat, 'econ'); 
sig = diag(S);

% % Create covariance matrix
% Cx = 1/(trimLen-1)*(firstMat*firstMat');
% [v,d] =  eig(Cx);
% % Cx = v * d * v'
% % We see that v contains principal components
% % First 
% theta = 90;
% R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% thirdcam = R'*[trimCam3x trimCam3y]';
% 
% firstMatUpdate = [trimCam1x'; trimCam1y'; trimCam2x'; trimCam2y'; thirdcam(1,:); thirdcam(2,:)];
% CxUpdate = 1/(trimLen-1)*(firstMatUpdate*firstMatUpdate');
% [vup,dup] = eig(CxUpdate);
% [Uup,Sup,Vup] = svd(firstMatUpdate, 'econ');
% sigup = diag(Sup);


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


% subplot(1,2,2)
% plot(sigup.^2/sum(sigup.^2),'ko','Linewidth',1)
% title('Sigma Energies of Ideal Case Rotated Matrix')
% xlabel('Sigma Index'); ylabel('Energy')
% figure
% hold on

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


prin_comp = U' * fourthMat;
plot(1:trimLen, prin_comp(3:5,:),'Linewidth',1)
xlabel('Time (frames)'); ylabel('Displacement (pixels)')
legend('3rd Principal Component', '4th Principal Component', ...
       '5th Principal Component', 'location', 'northeast')
title('Case 4: Three to Five Principal Component Projections')