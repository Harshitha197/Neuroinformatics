%% PhysioNet EEG Motor Imagery Dataset - Starter Pipeline
% Tested in MATLAB R2024a

% --- Setup ---
dataDir = 'physionet.org/files/eegmmidb/1.0.0/S001';
edfFile = fullfile(dataDir,'S001R03.edf');   % Example: Subject 1, Run 3

% --- Read EDF ---
[TT, ann] = edfread(edfFile);   % TT = timetable (signals), ann = timetable (annotations)
fs = 160;                       % sampling rate (Hz), fixed for this dataset

% --- Signals matrix ---
% Convert timetable variables (which might be cell) into numeric
if iscell(TT{:,1})
    % Each column is a cell containing a numeric vector
    nCh = width(TT);
    nSamp = numel(TT{1,1}{1});
    X = zeros(nCh, nSamp);
    for ch = 1:nCh
        X(ch,:) = TT{1,ch}{1};   % extract numeric vector from cell
    end
else
    % Already numeric
    X = table2array(TT).';
end

chNames = string(TT.Properties.VariableNames);

% --- Events ---
evOnsets = seconds(ann.Onset);          % in seconds
rawLabels = ann.Annotations;            % cell array of strings
expandedOnsets = [];
expandedLabels = [];

for k = 1:numel(rawLabels)
    labs = split(string(rawLabels{k}),',');
    for j = 1:numel(labs)
        if labs(j) ~= ""
            expandedOnsets(end+1,1) = evOnsets(k);
            expandedLabels{end+1,1} = strtrim(labs(j));
        end
    end
end

onsetSamples = round(expandedOnsets*fs)+1;

% --- Run number parsing ---
tok = regexp(edfFile,'R(\d+)\.edf','tokens','once');
runNum = str2double(tok{1});

% --- Label mapping (from dataset docs) ---
% Runs 3,4,7,8,11,12 → T1=Left Fist, T2=Right Fist
% Runs 5,6,9,10,13,14 → T1=Both Fists, T2=Both Feet
isLR = ismember(runNum,[3 4 7 8 11 12]);
isFB = ismember(runNum,[5 6 9 10 13 14]);

ixT1 = find(strcmp(expandedLabels,"T1"));
ixT2 = find(strcmp(expandedLabels,"T2"));
t1Samples = onsetSamples(ixT1);
t2Samples = onsetSamples(ixT2);

% --- Preprocessing: Bandpass + Notch + CAR ---
bp = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',1,'HalfPowerFrequency2',40,'SampleRate',fs);
notch = designfilt('bandstopiir','FilterOrder',4, ...
    'HalfPowerFrequency1',49,'HalfPowerFrequency2',51,'SampleRate',fs);

Xf = filtfilt(bp,X.').';
Xf = filtfilt(notch,Xf.').';
Xref = Xf - mean(Xf,1);   % Common average reference

% --- Epoching function ---
function epochs = makeEpochs(Xref, inds, swin)
    valid = inds(inds+swin(1)>0 & inds+swin(2)<=size(Xref,2));
    nValid = numel(valid);
    nCh = size(Xref,1);
    epochLen = swin(2)-swin(1)+1;
    epochs = zeros(nCh, epochLen, nValid);
    for k = 1:nValid
        s = valid(k);
        epochs(:,:,k) = Xref(:, s+swin(1):s+swin(2));
    end
end

% --- Build epochs around events ---
win = [-0.5 2.5];           % window in seconds
swin = round(win*fs);       % window in samples
epochLen = swin(2)-swin(1)+1;

epochsT1 = makeEpochs(Xref, t1Samples, swin);
epochsT2 = makeEpochs(Xref, t2Samples, swin);

% --- Baseline correction ---
preIdx = 1:(-swin(1));
epochsT1 = epochsT1 - mean(epochsT1(:,preIdx,:),2);
epochsT2 = epochsT2 - mean(epochsT2(:,preIdx,:),2);
want = ["C3","Cz","C4"];
idx = find(contains(chNames, want, 'IgnoreCase', true));
labels = chNames(idx);

ERP_T1 = squeeze(mean(epochsT1(idx,:,:),3));   % [channels × samples]
ERP_T2 = squeeze(mean(epochsT2(idx,:,:),3));   % [channels × samples]

nSamp = size(ERP_T1,2);
tEpoch = linspace(win(1)/fs, win(2)/fs, nSamp);
% {ERP.chanlocs.labels}

figure;
subplot(2,1,1);
plot(tEpoch, ERP_T1.'); xline(0,'k--');
title('ERP T1'); xlabel('Time (s)'); ylabel('µV'); 
legend(labels, "Interpreter","none");

subplot(2,1,2);
plot(tEpoch, ERP_T2.'); xline(0,'k--');
title('ERP T2'); xlabel('Time (s)'); ylabel('µV'); 
legend(labels, "Interpreter","none");


   % show first 10 channel labels

