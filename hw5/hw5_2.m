%% Compare Foreground, Background, Original
image=reshape(u_dmd_fg(:,105),540,960);
image=uint8(real(image));
imshow(image);drawnow
%% Background
image = reshape(u_dmd_bg(:,105),540,960);
image = uint8(real(image));
imshow(image); drawnow
%% Original
imshow(orig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Monte Carlo
v1 = VideoReader('monte_carlo_low.mp4');
% dt = 1; % dt is 1 frame
iter = 0;
col_images = zeros(518400,379);
while hasFrame(v1)
   iter = iter + 1;
   frame = readFrame(v1);
   frame = rgb2gray(frame);
   
   if iter == 50
       orig = frame;
   end
   
   sz = size(frame,1) * size(frame,2);
   image = reshape(frame(:,:),sz,1);
   col_images(:,iter) = image;
end
t = 1:iter;
dt = t(2) - t(1);
%% Get X1, X2, Perform SVD
X1 = col_images(:, 1:end-1);
X2 = col_images(:, 2:end);
[U,Sigma,V] = svd(X1,'econ');
sig = diag(Sigma);
%% Look at Singular Values
figure(1)
% subplot(1,2,1)
plot(sig(1:50).^2/sum(sig(1:50).^2),'ko','Linewidth',1)
xlabel('Sigma Index'); ylabel('Energy')
title('Sigma Energies')

figure(2)
% first 50 singular values
plot(sig(1:50), 'ko')
%%
mode = 20;
% use some subset of U V and Sigma to get low rank approx 
S = U(:,1:mode)'*X2*V(:,1:mode)*diag(1./diag(Sigma(1:mode,1:mode)));
[eV, D] = eig(S);
mu = diag(D);
Phi = U(:,1:mode)*eV; % dmd modes in cols
%% DMD Spectra
% plot omega in the complex plane - real vs imaginary of omega
%   pick some threshold to see which omega is close to magnitude 0
%   proceed with omegas that are part of background
omega = log(mu)/dt;

