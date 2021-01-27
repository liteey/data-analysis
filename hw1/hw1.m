% Clean workspace
clear; close all; clc

load subdata.mat

L = 10;
n = 64;
x2 = linspace(-L,L,n+1); x = x2(1:n); y = x; z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);

[X,Y,Z] = meshgrid(x,y,z);
[Kx,Ky,Kz] = meshgrid(ks,ks,ks);

% average in frequency domains
Utnave = zeros(n,n,n);
for j = 1:49
    Un(:,:,:) = reshape(subdata(:,j),n,n,n);
    Utn = fftn(Un);
    Utnave = Utnave + Utn;
    
end
Utnave = fftshift(Utnave)/49;
M = max(abs(Utnave),[],'all');
% Plot central frequency
isosurface(Kx,Ky,Kz,abs(Utnave)/M,0.7);
xlabel('Kx')
ylabel('Ky')
zlabel('Kz')
%%
% Create the Gaussian Filter
[Max, Ind] = max(abs(Utnave(:)));
[Ix, Iy, Iz] = ind2sub(size(Utnave), Ind);
Kx0 = Kx(Ix, Iy, Iz);  %  5.34
Ky0 = Ky(Ix, Iy, Iz);  % -6.91
Kz0 = Kz(Ix, Iy, Iz);  %  2.19
tau = 1.5;  % Can change
filter = exp(-tau * ((Kx - Kx0).^2 + (Ky - Ky0).^2 + (Kz - Kz0).^2));
isosurface(Kx,Ky,Kz,filter,0.7);
xlabel('Kx')
ylabel('Ky')
zlabel('Kz')
%%
% Find and plot path of the submarine
P = zeros(3,49);
for j = 1:49
    Un(:,:,:)= reshape(subdata(:,j),n,n,n);
    Unt = fftshift(fftn(Un));
    Untf = filter .* Unt;
    Unf = ifftn(fftshift(Untf));
    [Max, Ind] = max(abs(Unf(:)));
    [Ixx, Iyy, Izz] = ind2sub(size(abs(Unf)), Ind);
    P(:,j) = [Ixx, Iyy, Izz];
end
plot3(P(1,:), P(2,:), P(3,:));
xlabel('x')
ylabel('y')
zlabel('z')
