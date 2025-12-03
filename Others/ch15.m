clear; clc;
load('sampleEEGdata.mat'); 

if ~exist('data','var') || ~exist('srate','var')
    srate   = 1000;              
    t       = 0:1/srate:10;       
    foi     = 10;                
    nChans  = 32;
    data    = randn(nChans, numel(t))*0.2; 
    data(5,:) = data(5,:) + sin(2*pi*foi*t).*exp(-(t-5).^2/2); 
end

elecIdx = 5;           
foi     = 10;          
sig = data(elecIdx, :);
nTime = length(sig);
time  = (0:nTime-1)/srate;  
%  1) Complex Morlet wavelet convolution (power)
nCycles = 7; 
t_wav   = -1:1/srate:1;  

sigma_t = nCycles / (2*pi*foi);
cmw = exp(1i*2*pi*foi.*t_wav) .* exp(-t_wav.^2./(2*sigma_t^2));
cmw = cmw ./ sqrt(sum(abs(cmw).^2));

% Convolution (same length as signal)
convResult = conv(sig, cmw, 'same');

% Power time series = magnitude squared of analytic signal
power_wavelet = abs(convResult).^2;

%  2) Filter-Hilbert method (power)

fLow  = max(0.1, foi - 2);   
fHigh = foi + 2;             

filterOrder = 4;  

Wn = [fLow fHigh] / (srate/2);
[b, a] = butter(filterOrder, Wn, 'bandpass');

sig_filt = filtfilt(b, a, sig);
analytic_sig = hilbert(sig_filt);
power_hilbert = abs(analytic_sig).^2;

%  3) Short-time FFT / Spectrogram (power)
winLength_sec = 0.5;                           
winLength_smp = round(winLength_sec * srate);  
overlap_smp = round(0.5 * winLength_smp);

nFFT = 2^nextpow2(winLength_smp);
[S, F, T] = spectrogram(sig, winLength_smp, overlap_smp, nFFT, srate);

% Power spectrogram
P = abs(S).^2;  

[~, foi_idx] = min(abs(F - foi));
power_stft = P(foi_idx, :);  

figure;
set(gcf,'Position',[100 100 900 600])

% 1) Complex wavelet power
subplot(3,1,1);
plot(time, power_wavelet, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Power');
title(sprintf('Complex wavelet convolution power @ %g Hz (Electrode %d)', foi, elecIdx));
grid on;

% 2) Filter-Hilbert power
subplot(3,1,2);
plot(time, power_hilbert, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Power');
title(sprintf('Filter-Hilbert power @ %g Hz (Electrode %d)', foi, elecIdx));
grid on;

% 3) STFT power
subplot(3,1,3);
plot(T, power_stft, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Power');
title(sprintf('Short-time FFT power @ %g Hz (Electrode %d)', foi, elecIdx));
grid on;

sgtitle(sprintf('Power over time at %g Hz, electrode %d: comparison of 3 methods', foi, elecIdx));