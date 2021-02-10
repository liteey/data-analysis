%% Part 1 GNR Guitar Frequencies
clear; close all; clc;
[y1, Fs1] = audioread('GNR.m4a');
tr_gnr = length(y1)/Fs1; % record time in seconds

% Setting up Time and Frequency Domain
L = tr_gnr;
n = length((1:length(y1))/Fs1);
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% Turn y into a row vector
y = y1';

% Experimented with a to get a value that shows the frequencies well
a = 500; % 
% Go across time using tau
tau = 0:0.1:L;
filt_Yft_spec = zeros(length(y), length(tau));
for j = 1:length(tau)
    % Create the Gabor Transform at t = tau
    filter = exp(-a*(t - tau(j)).^2);
    % Apply the filter
    Yf = filter .* y;
    % Go into Frequency space
    Yft = fft(Yf);

    % Find the max value in k
    [Max, Ind] = max(abs(Yft));
    [Max_Ind] = ind2sub(size(Yft), Ind);
    Max_Val = abs(k(Max_Ind));
    
    % Create the Gaussian Filter, centered around the max value of k in
    % t = tau
    fft_tau = 0.0001;
    fft_filt = exp(-fft_tau*(k - Max_Val).^2);
    filt_Yft = fft_filt .* Yft;
   
    filt_Yft_spec(:,j) = (fftshift(abs(filt_Yft)));
end

pcolor(tau, ks, abs(filt_Yft_spec));
shading interp
title("Guns N' Roses - Sweet Child O' Mine Guitar Frequencies")
% Set up the ylimit (which is frequency) to be between 200 and 850 Hz
% Because that is where the observed frequencies are
set(gca, 'ylim', [200 850], 'fontsize', 11)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')

%% Part 2 - Floyd Bass Fequencies
clear; close all; clc;
% figure(2)
[y2, Fs2] = audioread('Floyd.m4a');
tr_floyd = length(y2)/Fs2;

% Create time and frequency domain
L = tr_floyd;
n = length((1:length(y2)-1)/Fs2);
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% y's length is initially odd and because we had to divide n by 2 in
% frequency domain, we have to lose the last column of data so the vectors
% are the same length.
y = y2(1:length(y2)-1)';

% Experimented with a to get a value that shows the frequencies well
a = 100;
% Go across time using tau
tau = 0:0.5:10;
filt_Yft_spec = zeros(length(y), length(tau));
floyd_ind = zeros(1,length(tau));
floyd_freq = zeros(1,length(tau));
for j = 1:length(tau)
    % Create the Gabor Filter
    gabor = exp(-a*(t - tau(j)).^2);
    Yf = gabor .* y;
    Yft = fft(Yf);
    
    % Find the maximum
    [Max, Ind] = max(abs(Yft));
    [Max_Ind] = ind2sub(size(Yft), Ind);
    Max_Val = abs(k(Max_Ind));
    
    % Create the Gaussian Filter with constant 0.1 and filter around the
    % maximum that we found earlier.
    fft_tau = 0.1;
    fft_filt = exp(-fft_tau*(k - Max_Val).^2);
    filt_Yft = fft_filt .* Yft;

    filt_Yft_spec(:,j) = (fftshift(abs(filt_Yft)));

end

pcolor(tau,ks,log(abs(filt_Yft_spec)+1));
shading interp
title("Pink Floyd - Comfortably Numb Bass Frequencies")
% Set up the ylimit (which is frequency) to be between 65 and 150 Hz
% Because we notice that is where all the bass lies in.
set(gca, 'ylim', [65 150], 'fontsize', 11)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')
%% Part 3 - Floyd Guitar Frequencies
% clear; close all; clc;
% figure(2)
[y2, Fs2] = audioread('Floyd.m4a');
tr_floyd = length(y2)/Fs2;

L = tr_floyd;
n = length((1:length(y2)-1)/Fs2);
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% y's length is initially odd and because we had to divide n by 2 in
% frequency domain, we have to lose the last column of data so the vectors
% are the same length.
y = y2(1:length(y2)-1)';

pre_filt_Yft = fftshift(fft(y));
pre_filt_Yft(abs(pre_filt_Yft) >= 300) = 0;
pre_filt_Y = ifft(ifftshift(pre_filt_Yft));
    
% Experimented with a to get a value that shows the frequencies well
a = 100;
% Go across time using tau
tau = 30:0.5:L;
filt_Yft_spec = zeros(length(y), length(tau));
for j = 1:length(tau)
    % Create the Gabor Filter
    gabor = exp(-a*(t - tau(j)).^2);
    Yf = gabor .* pre_filt_Y;
    Yft = fft(Yf);
    
    % Find the maximum
    [Max, Ind] = max(abs(Yft));
    [Max_Ind] = ind2sub(size(Yft), Ind);
    Max_Val = abs(k(Max_Ind));
    
    % Create the Gaussian Filter with constant 0.1 and filter around the
    % maximum that we found earlier.
    fft_tau = 0.1;
    fft_filt = exp(-fft_tau*(k - Max_Val).^2);
    filt_Yft = fft_filt .* Yft;

    filt_Yft_spec(:,j) = (fftshift(abs(filt_Yft)));

end

pcolor(tau,ks,filt_Yft_spec);
shading interp
title("Pink Floyd - Comfortably Numb Presumed Guitar Frequencies")
% Set up the ylimit (which is frequency) to be between 0 and 150 Hz
% Because that is the range limit of a guitar
set(gca, 'ylim', [50 1000], 'fontsize', 11)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')