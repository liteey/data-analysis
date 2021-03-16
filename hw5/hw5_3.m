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