clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do not edit this part
%geneate signal to analyse
fs = 100;                                % sample frequency (Hz)
t = 0:1/fs:10-1/fs;                      % 10 second span time vector
x = (1.3)*sin(2*pi*15*t) ...             % 15 Hz component
  + (1.7)*sin(2*pi*40*(t-2)) ...         % 40 Hz component
  + 2.5*gallery('normaldata',size(t),4); % Gaussian noise;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(t);
X = fft(x); % FFT computation
X = fftshift(x);
X = X/N;

dt = mean(diff(t)); % sample spacing
df = 1/(N*dt); % frequency spacing
fi = (0:(N-1)) - floor(N/2); % unfolded index
f = df * fi; % frequency vector

power = abs(X) .^ 2;
plot(f, power); % power vs frequency
title('Power vs Frequency');
xlabel('Frequency');
ylabel('Power');