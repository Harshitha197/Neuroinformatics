clc
clear
close all

% Load
load('sampleEEGdata.mat');
signal = EEG.data(10,:,5); % 10th chan data for the 5th trial 
t = (0:length(signal)-1) / EEG.srate;

% Plot
figure;
plot(t, signal);
xlim([0 max(t)])
xlabel('Time (s)');
ylabel('Amplitude (µV)');
title('Raw EEG Signal (Single Channel)');

% Do a DFT
N = length(signal);
n = 0:N-1;
k = n';
W = exp(-1i*2*pi/N * (k * n));

X_manual = W * signal';  % DFT result
frequencies = (0:N-1) * (EEG.srate/N); 

% Do a FFT
X_fft = fft(signal)';
% X_fft = fft(signal,5*length(frequencies))'; % Hmm, can I really trust the values?
disp(length(frequencies));
disp(length(abs(X_fft)))

figure;
plot(frequencies, abs(X_manual));
hold on
plot(frequencies, X_fft);
xlim([0 frequencies(round(end))])
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('DFT', 'FFT');
title('FFT of EEG Signal (Custom fn)');

% Gaussian Kernel
time = -1:1/EEG.srate:1;
s = 5/(2*pi*30);
kernel = exp((-time.^2)/(2*s^2))/30;
kernel = kernel / sum(kernel);

figure;
plot(kernel, '-', 'LineWidth', 2);
xlabel('Time Points'); ylabel('Amplitude');
title('Kernel: Gaussian');

% Manual
conv_length = length(signal) + length(kernel) - 1;
conv_result = zeros(1, length(signal) + length(kernel) - 1);
sig_padded = [signal, zeros(1, length(kernel)-1)];
kernel_padded = [kernel, zeros(1, length(signal)-1)];

for n = 1:conv_length
    for k = 1:length(kernel)
        if (n-k+1 > 0) && (n-k+1 <= length(signal))
            conv_result(n) = conv_result(n) + kernel(k) * signal(n-k+1);
        end
    end
end

% DTFT
N_dtft = conv_length;
omega = 2 * pi * (0:N_dtft-1) / N_dtft;

X_dtft = zeros(1, N_dtft);
K_dtft = zeros(1, N_dtft);

for w_idx = 1:N_dtft
    w = omega(w_idx);
    X_dtft(w_idx) = sum(signal .* exp(-1i * w * (0:length(signal)-1)));
    K_dtft(w_idx) = sum(kernel .* exp(-1i * w * (0:length(kernel)-1)));
end

Y_dtft = X_dtft .* K_dtft;

conv_dtft = zeros(1, N_dtft);
for n = 1:N_dtft
    conv_dtft(n) = sum(Y_dtft .* exp(1i * omega * (n-1))) / N_dtft;
end
conv_dtft = real(conv_dtft);

% FFT and IFFT
N_fft = conv_length;
fft_sig = fft([signal, zeros(1, N_fft - length(signal))], N_fft);
fft_kernel = fft([kernel, zeros(1, N_fft - length(kernel))], N_fft);
conv_fft = ifft(fft_sig .* fft_kernel);

% Trim
conv_result_same = conv_result(floor(length(kernel)/2)+1:end-floor(length(kernel)/2));
conv_fft_same = real(conv_fft(floor(length(kernel)/2)+1:end-floor(length(kernel)/2)));
conv_dtft_same = real(conv_dtft(floor(length(kernel)/2)+1:end-floor(length(kernel)/2)));

% Plots
figure;
plot(signal, '-', 'LineWidth', 1);
hold on;
plot(conv_result_same, '-', 'LineWidth', 1);
hold on;
plot(real(conv_fft_same), '-', 'LineWidth', 1.5);
hold on;
plot(real(conv_dtft_same), '-', 'LineWidth', 1.5);
xlabel('Time Points'); ylabel('Amplitude');
title('Convolution Result (Manual vs FFT)');
legend('Raw EEG', 'Manual Convolution', 'FFT Convolution', 'DTFT Convolution'); % All is same ...

figure;
N_fft = length(conv_fft);
frequencies_conv = (0:N_fft-1) * (EEG.srate / N_fft);
plot(frequencies_conv, abs(fft_sig .* fft_kernel)); % Convolution in time domain is multiplication of their transforms in freq domain
hold on
plot(frequencies, abs(X_fft))
xlim([0 frequencies_conv(round(end/2))])
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT of Convolved EEG Signal'); % See how Guassian is a low pass filter

% Can you think of a kernel that will double up as a high pass filter?

% try doing Ex. 3 of p. 19 of Cohen's book.