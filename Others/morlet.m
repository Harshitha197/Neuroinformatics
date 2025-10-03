clear
close all
clc

m = 100; % Wavelet number; increase this to gain frequency precision at cost of temporal ...
% Increase to perhaps 20, in the next run ...
% Then, down to 2?

fs = 1000; % sampling frequency (Hz)
t = 0:1/fs:100;
signal = zeros(size(t));

signal(t <= 30) = sin(2*pi*2*t(t <= 30)); % 0–30 s → 2 Hz sine
signal(t > 30 & t <= 75) = sin(2*pi*4*t(t > 30 & t <= 75)); % 30–75 s → 4 Hz sine
signal(t > 75) = sin(2*pi*6*t(t > 75)); % 75–100 s → 6 Hz sine

figure;
plot(t, signal)
title('Signal');
w = morlet_wavelet(10, 6, t-mean(t));
hold on
plot(t, real(w));
title('Signal + Morlet wavelet (real; C = 10 Hz)');

freqs = 1:1:40;
power = zeros(length(freqs), length(t));

for fi = 1:length(freqs)
    f = freqs(fi);
    sigma = m/(2*pi*f);
    tw = -6*sigma : 1/fs : 6*sigma;   
    % time vector for morlet; scale s.t. the entire wavelet can be accomodated 
    w = morlet_wavelet(f, m, tw);
    conv_res = conv(signal, w, 'same');
    power(fi,:) = abs(conv_res).^2;
end

figure;
imagesc(t, freqs, power);
axis xy; xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Time-Frequency Power with COI'); colorbar;

function w = morlet_wavelet(f0, m, t)
    sigma = m / (2*pi*f0);
    w = exp(2*1i*pi*f0*t) .* exp(-t.^2/(2*sigma^2));
    w = w ./ sqrt(sum(abs(w).^2)) / sqrt(sigma); % normalize energy and scale
end