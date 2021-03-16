clear; close all; clc;
%%
v1 = VideoReader('ski_drop_low.mp4');
% dt = 1; % dt is 1 frame
iter = 0;
col_images = zeros(518400,454);
while hasFrame(v1)
   iter = iter + 1;
   frame = readFrame(v1);
   frame = rgb2gray(frame);
   
   if iter == 105
       orig = frame;
   end
   
   sz = size(frame,1) * size(frame,2);
   image = reshape(frame(:,:),sz,1);
   col_images(:,iter) = image;
end
t = 1:iter;
dt = t(2) - t(1);
%% get x1 and x2
X1 = col_images(:, 1:end-1);
X2 = col_images(:, 2:end);
%% do svd, get eigenvals and eigenvectors
[U,Sigma,V] = svd(X1,'econ');
sig = diag(Sigma);
%%
energy = 0;
modes = 0;
while energy < 0.9
    modes = modes + 1;
    energy = energy + (sig(modes).^2/sum(sig.^2));

    
end
% modes = modes - 1;

%%
figure(1)
% subplot(1,2,1)
plot(sig.^2/sum(sig.^2),'ko','Linewidth',1)
xlabel('Sigma Index'); ylabel('Energy')
title('Sigma Energies')
%%
mode = 20;
% use some subset of U V and Sigma to get low rank approx 
S = U(:,1:mode)'*X2*V(:,1:mode)*diag(1./diag(Sigma(1:mode,1:mode)));
[eV, D] = eig(S);
mu = diag(D);
Phi = U(:,1:mode)*eV; % dmd modes in cols
%%
% first 50 singular values
plot(sig(1:50), 'ko')

%% DMD Spectra
% plot omega in the complex plane - real vs imaginary of omega
%   pick some threshold to see which omega is close to magnitude 0
%   proceed with omegas that are part of background
omega = log(mu)/dt;

figure
hold on
plot([-0.6 0.3], [0 0], 'k')
plot([0 0], [-0.1 0.1],'k')
plot(real(omega), imag(omega),'o');
title('Real vs Complex Omega Values')
% we see that some points close to real and others close to imaginary
%%
% plot of singular values with we find svd
% plot of dmd eigenvalues in complex plane
% final results: how well of separating foreground and background:
%   show some sample frames, show full frame, show bg, and show bg
thresh = 0.001;
bg = find(abs(omega) < thresh);
omega_bg = omega(bg);
phi_bg = Phi(:,bg);

%% Compute DMD construction of background
y0 = phi_bg\X1(:,1);


% y0 = Phi\X1(:,1);
u_modes = zeros(length(y0),iter);
%%
for j = 1:iter
   u_modes(:,j) = y0.*exp(omega_bg*t(j)); 
end
u_dmd_bg = phi_bg*u_modes;
%% Subtract background to get foreground
u_dmd_fg = col_images - abs(u_dmd_bg);
ind = find(u_dmd_fg < 0);
X_bgr = u_dmd_bg;
X_bgr(ind) = u_dmd_bg(ind) + u_dmd_fg(ind);
X_fgr = u_dmd_fg;
X_fgr(ind) = 0;
%% Show Foreground



% for j=1:iter
%    image=reshape(u_dmd_fg(:,j),540,960);
%    image=uint8(real(image));
%    imshow(image);drawnow
% end


%% Show Background
remake = zeros(540,960,1,iter);
for j = 1:iter
    image = reshape(u_dmd_bg(:,j),540,960);
    image = uint8(real(image));
    imshow(image); drawnow
%     remake(:,:,1,j) = image;
end
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

figure
hold on
plot([-0.3 0.1], [0 0], 'k')
plot([0 0], [-0.4 0.4],'k')
plot(real(omega), imag(omega),'o');
title('Real vs Complex Omega Values')
% we see that there are more points that are off the axes unlike skidrop
%%
thresh = 0.001;
bg = find(abs(omega) < thresh);
omega_bg = omega(bg);
phi_bg = Phi(:,bg);
%% Compute DMD construction of background
y0 = phi_bg\X1(:,1);
u_modes = zeros(length(y0),iter);
for j = 1:iter
   u_modes(:,j) = y0.*exp(omega_bg*t(j)); 
end
u_dmd_bg = phi_bg*u_modes;
%% Subtract background to get foreground
u_dmd_fg = col_images - abs(u_dmd_bg);
ind = find(u_dmd_fg < 0);
X_bgr = u_dmd_bg;
X_bgr(ind) = u_dmd_bg(ind) + u_dmd_fg(ind);
X_fgr = u_dmd_fg;
X_fgr(ind) = 0;
%% Show foreground
for j=1:iter
   image=reshape(u_dmd_fg(:,j),540,960);
   image=uint8(real(image));
   imshow(image);drawnow
end
%% Compare fg, bg and original
% fg
image=reshape(u_dmd_fg(:,50),540,960);
image=uint8(real(image));
imshow(image);drawnow
%% bg
image = reshape(u_dmd_bg(:,50),540,960);
image = uint8(real(image));
imshow(image); drawnow
%% Original
imshow(orig)