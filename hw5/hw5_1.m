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
