
clear; clc;
load sampleEEGdata.mat   

data = EEG.data;
times = EEG.times;         
chanlocs = EEG.chanlocs;   
srate = EEG.srate;         
nChans = EEG.nbchan;
nTimes = length(times);

freqsOfInterest = [6 10 20];        
timePointsToPlot = [200 400 600 800 1000];   
baselineWindow = [-500 -200];         
nCycles = 6;                       

nFreqs = length(freqsOfInterest);

%   STEP 1: Wavelet time-frequency POWER: chans × freqs × times
power_tf = zeros(nChans, nFreqs, nTimes);

for fi = 1:nFreqs
    freq = freqsOfInterest(fi);

    % Create complex Morlet wavelet
    t_wav = -1:1/srate:1;
    s = nCycles/(2*pi*freq);
    wavelet = exp(1i*2*pi*freq.*t_wav) .* exp(-t_wav.^2/(2*s^2));
    wavelet = wavelet ./ sqrt(sum(abs(wavelet).^2));

    % Convolve each channel
    for ch = 1:nChans
        sig = double(squeeze(data(ch,:,1))); 
        convRes = conv(sig, wavelet, 'same');
        power_tf(ch,fi,:) = abs(convRes).^2;
    end
end

%   STEP 2: Baseline normalization (dB)
[~, b1] = min(abs(times - baselineWindow(1)));
[~, b2] = min(abs(times - baselineWindow(2)));

baselinePower = mean(power_tf(:,:,b1:b2), 3);  
power_dB = zeros(size(power_tf));

for fi = 1:nFreqs
    power_dB(:,fi,:) = 10 * log10( bsxfun(@rdivide, squeeze(power_tf(:,fi,:)), baselinePower(:,fi)) );
end

%   STEP 3: Create topographical maps
%   5 time points × 2 rows per frequency
for fi = 1:nFreqs
    
    freq = freqsOfInterest(fi);
    pow_f = squeeze(power_tf(:,fi,:));   
    pow_db = squeeze(power_dB(:,fi,:));  

    clim_raw = [min(pow_f(:)) max(pow_f(:))];
    clim_db  = [min(pow_db(:)) max(pow_db(:))];

    figure('Position',[100 50 1300 550]);
    sgtitle(sprintf('Topographical Power @ %d Hz', freq));

    for ti = 1:length(timePointsToPlot)
        [~, tidx] = min(abs(times - timePointsToPlot(ti)));

        %% ---- ROW 1: RAW POWER ----
        subplot(2, length(timePointsToPlot), ti)
        topoplot(pow_f(:,tidx), chanlocs, 'maplimits', clim_raw);
        title(sprintf('%d ms (raw)', times(tidx)));
        colorbar

        %% ---- ROW 2: BASELINE NORMALIZED POWER ----
        subplot(2, length(timePointsToPlot), length(timePointsToPlot) + ti)
        topoplot(pow_db(:,tidx), chanlocs, 'maplimits', clim_db);
        title(sprintf('%d ms (dB)', times(tidx)));
        colorbar
    end
end