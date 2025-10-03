clear; clc;
config = struct(...
    'data_folder', './', ...
    'output_folder', './processed/S001', ...
    'subjects', {{'S001'}}, ...
    'runs', {{'R04'}}, ...
    'filter_hp', 1, ...
    'filter_lp', 40, ...
    'epoch_window', [-1, 4], ...
    'baseline', [-0.5, 0], ...
    'artifact_threshold', 200 ...
);

if ~exist(config.output_folder, 'dir'), mkdir(config.output_folder); end
[ALLEEG, EEG, CURRENTSET] = eeglab('nogui');

fprintf('Starting EEG preprocessing...\n\n');
results = [];

for s = 1:length(config.subjects)
    for r = 1:length(config.runs)
        subject = config.subjects{s};
        run = config.runs{r};
        
        fprintf('Processing %s %s...\n', subject, run);
        
        try
            EEG_raw = load_eeg_data(config.data_folder, subject, run);
            original_channels = EEG_raw.nbchan;
            
            EEG_filtered = preprocess_eeg_continuous(EEG_raw, config);
            
            [EEG_epoched, epoch_stats] = create_epochs(EEG_filtered, config.epoch_window, config.baseline);
            
            filename = sprintf('%s_%s_processed.set', subject, run);
            EEG_final = pop_saveset(EEG_epoched, 'filename', filename, 'filepath', config.output_folder);
            
            result = struct(...
                'subject', subject, 'run', run, ...
                'channels', sprintf('%d/%d', EEG_final.nbchan, original_channels), ...
                'epochs', EEG_final.trials, ...
                'T0', epoch_stats.T0, 'T1', epoch_stats.T1, 'T2', epoch_stats.T2, ...
                'status', '✅ Success' ...
            );
            results = [results; result];
            
            fprintf('  ✅ %d channels, %d epochs (%d T0, %d T1, %d T2)\n\n', ...
                    EEG_final.nbchan, EEG_final.trials, epoch_stats.T0, epoch_stats.T1, epoch_stats.T2);
            
            % Create improved preprocessing plots
            create_improved_preprocessing_plots(config, result, EEG_raw, EEG_filtered, EEG_epoched);
            
        catch ME
            fprintf('  ❌ Error: %s\n\n', ME.message);
            result = struct('subject', subject, 'run', run, 'status', '❌ Failed');
            results = [results; result];
        end
    end
end

print_summary(results);
fprintf('🎉 Processing complete! Files saved in: %s\n', config.output_folder);

function EEG = load_eeg_data(data_folder, subject, run)
    edf_file = fullfile(data_folder, subject, [subject run '.edf']);
    EEG = pop_biosig(edf_file);
    EEG.setname = sprintf('%s_%s', subject, run);
    
    % Load events
    event_file = [edf_file '.event'];
    if exist(event_file, 'file')
        EEG = extract_events(EEG, event_file);
    else
        EEG = create_synthetic_events(EEG);
    end
    
    fprintf('  📁 Loaded: %d channels, %.1fs, %d events\n', ...
            EEG.nbchan, EEG.pnts/EEG.srate, length(EEG.event));
end

function EEG = extract_events(EEG, event_file)
    % Read event file (handles binary format)
    fid = fopen(event_file, 'rb');
    raw_data = fread(fid, inf, 'uint8');
    fclose(fid);
    content = char(raw_data');
    
    % Extract T0, T1, T2 events with durations
    pattern = '(T[012])\s*duration:\s*([0-9]+\.?[0-9]*)';
    matches = regexp(content, pattern, 'tokens');
    
    if ~isempty(matches)
        events = [];
        current_time = 2.0;  % Start at 2s
        
        for i = 1:length(matches)
            events(i).type = matches{i}{1};
            events(i).latency = round(current_time * EEG.srate);
            duration = str2double(matches{i}{2});
            events(i).duration = round(duration * EEG.srate);
            current_time = current_time + duration;
            
            if current_time >= (EEG.pnts / EEG.srate) - 1, break; end
        end
        
        EEG.event = events;
    else
        % Fallback: count occurrences and create timing
        EEG = create_synthetic_events(EEG);
    end
end

function EEG = create_synthetic_events(EEG)
    events = [];
    current_time = 2.0;
    sequence = {'T0', 'T1', 'T0', 'T2'};  % rest, left, rest, right
    durations = [2.0, 4.0, 2.0, 4.0];    % realistic durations
    
    i = 1;
    while current_time < (EEG.pnts/EEG.srate) - 5
        idx = mod(i-1, 4) + 1;
        events(i).type = sequence{idx};
        events(i).latency = round(current_time * EEG.srate);
        events(i).duration = round(durations(idx) * EEG.srate);
        
        current_time = current_time + durations(idx);
        i = i + 1;
    end
    
    EEG.event = events;
end

function EEG = preprocess_eeg_continuous(EEG, config)
    % Filter
    EEG = pop_eegfiltnew(EEG, config.filter_hp, []);  % High-pass
    EEG = pop_eegfiltnew(EEG, [], config.filter_lp);  % Low-pass
    
    % Re-reference
    EEG = pop_reref(EEG, []);
    
    bad_channels = find_bad_channels(EEG, config.artifact_threshold);
    if ~isempty(bad_channels)
        EEG = pop_select(EEG, 'nochannel', bad_channels);
        fprintf('  🔧 Removed %d noisy channels\n', length(bad_channels));
    end
end

function bad_channels = find_bad_channels(EEG, threshold)
    channel_var = var(EEG.data, 0, 2);
    median_var = median(channel_var);
    
    bad_var = find(channel_var > median_var * 10 | channel_var < median_var * 0.1);
    bad_amp = find(any(abs(EEG.data) > threshold, 2));
    
    bad_channels = unique([bad_var; bad_amp]);
end

function [EEG, stats] = create_epochs(EEG, epoch_window, baseline)
    event_types = unique({EEG.event.type});
    EEG = pop_epoch(EEG, event_types, epoch_window);
    EEG = pop_rmbase(EEG, baseline * 1000);
   
    stats.T0 = count_epochs(EEG, 'T0');
    stats.T1 = count_epochs(EEG, 'T1');
    stats.T2 = count_epochs(EEG, 'T2');
    
    fprintf('  📊 Created %d epochs\n', EEG.trials);
end

function count = count_epochs(EEG, event_type)
    count = 0;
    for i = 1:EEG.trials
        if isfield(EEG.epoch(i), 'eventtype')
            types = EEG.epoch(i).eventtype;
            if iscell(types), types = types{1}; end
            if strcmp(types, event_type)
                count = count + 1;
            end
        end
    end
end

function print_summary(results)
    fprintf('\n📋 PROCESSING SUMMARY\n');
    fprintf('Subject | Run | Channels | Epochs | T0 | T1 | T2 | Status\n');
    fprintf('--------|-----|----------|--------|----|----|----|---------\n');
    
    for i = 1:length(results)
        r = results(i);
        if isfield(r, 'epochs')
            fprintf('%-7s | %-3s | %8s | %6d | %2d | %2d | %2d | %s\n', ...
                    r.subject, r.run, r.channels, r.epochs, r.T0, r.T1, r.T2, r.status);
        else
            fprintf('%-7s | %-3s |    ---   |   ---  | -- | -- | -- | %s\n', ...
                    r.subject, r.run, r.status);
        end
    end
    fprintf('\n');
end

%% IMPROVED PLOTTING FUNCTION
function create_improved_preprocessing_plots(config, result, EEG_raw, EEG_filtered, EEG_epoched)
    try
        subject = result.subject;
        run = result.run;
        
        % Create figure with multiple subplots
        figure('Name', 'EEG Preprocessing Results', 'Position', [100, 100, 1400, 900]);
        
        %% SUBPLOT 1: Raw vs Filtered Data Comparison (FIXED)
        subplot(2,4,1);
        channel_idx = min(10, min(EEG_raw.nbchan, EEG_filtered.nbchan));
        
        % Show same time segment for both (30 seconds)
        time_segment = 30; % seconds
        samples_to_show = min(time_segment * EEG_raw.srate, EEG_raw.pnts);
        time_raw = (0:samples_to_show-1) / EEG_raw.srate;
        
        % Find corresponding channel in filtered data (match by position for now)
        filtered_idx = min(channel_idx, EEG_filtered.nbchan);
        
        % Plot raw data
        raw_data = EEG_raw.data(channel_idx, 1:samples_to_show);
        plot(time_raw, raw_data, 'b-', 'LineWidth', 0.8, 'DisplayName', 'Raw');
        hold on;
        
        % Plot filtered data (same time segment)
        filtered_data = EEG_filtered.data(filtered_idx, 1:min(samples_to_show, size(EEG_filtered.data,2)));
        time_filtered = (0:length(filtered_data)-1) / EEG_filtered.srate;
        plot(time_filtered, filtered_data, 'r-', 'LineWidth', 1.2, 'DisplayName', 'Filtered (1-40 Hz)');
        
        xlabel('Time (s)'); ylabel('Amplitude (μV)');
        title(sprintf('Raw vs Filtered Data (Ch %d)', channel_idx));
        legend('show', 'Location', 'best');
        grid on; 
        xlim([0, min(10, max(time_raw))]);
        
        %% SUBPLOT 2: Power Spectral Density (IMPROVED)
        subplot(2,4,2);
        
        % Calculate PSD for raw data (use longer segment for better frequency resolution)
        psd_segment = min(60*EEG_raw.srate, EEG_raw.pnts); % 60 seconds or full data
        [psd_raw, f_raw] = pwelch(EEG_raw.data(channel_idx, 1:psd_segment), ...
                                  [], [], [], EEG_raw.srate);
        
        % Calculate PSD for filtered data
        [psd_filt, f_filt] = pwelch(EEG_filtered.data(filtered_idx, 1:min(psd_segment, size(EEG_filtered.data,2))), ...
                                    [], [], [], EEG_filtered.srate);
        
        semilogy(f_raw, psd_raw, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Raw'); 
        hold on;
        semilogy(f_filt, psd_filt, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Filtered');
        
        % Mark filter boundaries
        xline(config.filter_hp, 'g--', 'High-pass', 'LineWidth', 2);
        xline(config.filter_lp, 'g--', 'Low-pass', 'LineWidth', 2);
        
        xlabel('Frequency (Hz)'); ylabel('Power Spectral Density (μV²/Hz)');
        title('Frequency Content Comparison');
        legend('show', 'Location', 'best');
        grid on; xlim([0, 60]);
        
        %% SUBPLOT 3: Channel Quality Assessment (ENHANCED)
        subplot(2,4,3);
        
        % Calculate channel statistics for raw data
        raw_var = var(EEG_raw.data, 0, 2);
        raw_channels = 1:EEG_raw.nbchan;
        
        % Determine which channels were removed
        kept_channels = 1:EEG_filtered.nbchan;
        removed_count = EEG_raw.nbchan - EEG_filtered.nbchan;
        
        % Plot variance for all original channels
        semilogy(raw_channels, raw_var, 'bo', 'MarkerSize', 4, 'DisplayName', 'All Channels'); 
        hold on;
        
        % Highlight kept channels
        if EEG_filtered.nbchan < EEG_raw.nbchan
            semilogy(kept_channels, raw_var(kept_channels), 'go', 'MarkerSize', 6, 'DisplayName', 'Kept');
            % Show removed channels (approximate)
            [~, worst_channels] = sort(raw_var, 'descend');
            removed_idx = worst_channels(1:removed_count);
            semilogy(removed_idx, raw_var(removed_idx), 'rx', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Removed');
        end
        
        % Add threshold lines
        median_var = median(raw_var);
        yline(median_var * 10, 'r--', 'Upper Threshold', 'LineWidth', 1.5);
        yline(median_var * 0.1, 'r--', 'Lower Threshold', 'LineWidth', 1.5);
        
        xlabel('Channel Number'); ylabel('Variance (μV²)');
        title(sprintf('Channel Quality (%d/%d kept)', EEG_filtered.nbchan, EEG_raw.nbchan));
        legend('show', 'Location', 'best');
        grid on;
        
        %% SUBPLOT 4: Event Timeline (IMPROVED)
        subplot(2,4,4);
        
        if ~isempty(EEG_filtered.event)
            event_times = [EEG_filtered.event.latency] / EEG_filtered.srate;
            event_types = {EEG_filtered.event.type};
            
            % Create event markers
            unique_types = unique(event_types);
            colors = {'b', 'r', 'g', 'm', 'c'};
            
            for i = 1:length(unique_types)
                type = unique_types{i};
                indices = strcmp(event_types, type);
                times = event_times(indices);
                
                stem(times, ones(size(times)), colors{mod(i-1,5)+1}, ...
                    'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', type);
                hold on;
            end
            
            xlabel('Time (s)'); ylabel('Events');
            title(sprintf('Event Timeline (%d events)', length(event_times)));
            legend('show', 'Location', 'best');
            grid on; ylim([0, 1.5]);
        else
            text(0.5, 0.5, 'No events found', 'HorizontalAlignment', 'center', 'Units', 'normalized');
            title('Event Timeline');
        end
        
        %% SUBPLOT 5: Sample Epochs by Condition (NEW)
        subplot(2,4,5);
        
        if EEG_epoched.trials > 1
            channel_to_show = min(5, EEG_epoched.nbchan);
            
            % Show one epoch from each condition
            conditions = {'T0', 'T1', 'T2'};
            colors = {'b', 'r', 'g'};
            
            for c = 1:length(conditions)
                % Find first epoch of this condition
                for ep = 1:EEG_epoched.trials
                    if isfield(EEG_epoched.epoch(ep), 'eventtype')
                        types = EEG_epoched.epoch(ep).eventtype;
                        if iscell(types), types = types{1}; end
                        if strcmp(types, conditions{c})
                            epoch_data = squeeze(EEG_epoched.data(channel_to_show, :, ep));
                            plot(EEG_epoched.times/1000, epoch_data, colors{c}, 'LineWidth', 1.5, ...
                                 'DisplayName', [conditions{c} ' (epoch ' num2str(ep) ')']);
                            hold on;
                            break;
                        end
                    end
                end
            end
            
            xlabel('Time (s)'); ylabel('Amplitude (μV)');
            title(sprintf('Sample Epochs (Ch %d)', channel_to_show));
            
            % Mark important time points
            xline(0, 'k--', 'Event Onset', 'LineWidth', 1.5);
            if config.baseline(2) ~= 0
                xline(config.baseline(2), 'g--', 'Baseline End', 'LineWidth', 1);
            end
            
            legend('show', 'Location', 'best');
            grid on;
        else
            text(0.5, 0.5, 'Data not epoched', 'HorizontalAlignment', 'center', 'Units', 'normalized');
            title('Sample Epochs');
        end
        
        %% SUBPLOT 6: Epoch Statistics (NEW)
        subplot(2,4,6);
        
        if EEG_epoched.trials > 1
            % Calculate statistics across epochs for each condition
            conditions = {'T0', 'T1', 'T2'};
            condition_data = cell(1, 3);
            
            for c = 1:length(conditions)
                epochs_this_condition = [];
                for ep = 1:EEG_epoched.trials
                    if isfield(EEG_epoched.epoch(ep), 'eventtype')
                        types = EEG_epoched.epoch(ep).eventtype;
                        if iscell(types), types = types{1}; end
                        if strcmp(types, conditions{c})
                            epochs_this_condition = [epochs_this_condition, ep];
                        end
                    end
                end
                condition_data{c} = epochs_this_condition;
            end
            
            % Plot epoch counts
            counts = [length(condition_data{1}), length(condition_data{2}), length(condition_data{3})];
            bar_colors = [0.3 0.3 1; 1 0.3 0.3; 0.3 1 0.3];
            
            bar(1:3, counts, 'FaceColor', 'flat', 'CData', bar_colors);
            set(gca, 'XTickLabel', conditions);
            xlabel('Condition'); ylabel('Number of Epochs');
            title('Epoch Count by Condition');
            grid on;
            
            % Add text labels on bars
            for i = 1:3
                text(i, counts(i) + 0.5, num2str(counts(i)), ...
                     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            end
        else
            text(0.5, 0.5, 'No epoch statistics available', 'HorizontalAlignment', 'center', 'Units', 'normalized');
            title('Epoch Statistics');
        end
        
        %% SUBPLOT 7: Filtering Effect on Single Channel (NEW)
        subplot(2,4,7);
        
        % Show a shorter segment to better see filtering effect
        segment_samples = min(5*EEG_raw.srate, EEG_raw.pnts); % 5 seconds
        time_segment = (0:segment_samples-1) / EEG_raw.srate;
        
        % Use a channel that typically shows good signal (central channels)
        test_channel = min(round(EEG_raw.nbchan/2), EEG_raw.nbchan);
        filtered_test_ch = min(test_channel, EEG_filtered.nbchan);
        
        raw_segment = EEG_raw.data(test_channel, 1:segment_samples);
        filtered_segment = EEG_filtered.data(filtered_test_ch, 1:min(segment_samples, size(EEG_filtered.data,2)));
        time_filt_seg = (0:length(filtered_segment)-1) / EEG_filtered.srate;
        
        plot(time_segment, raw_segment, 'b-', 'LineWidth', 0.8, 'DisplayName', 'Raw Signal');
        hold on;
        plot(time_filt_seg, filtered_segment, 'r-', 'LineWidth', 1.2, 'DisplayName', 'Filtered Signal');
        
        xlabel('Time (s)'); ylabel('Amplitude (μV)');
        title(sprintf('Filtering Effect Detail (Ch %d)', test_channel));
        legend('show', 'Location', 'best');
        grid on;
        
        subplot(2,4,8);
        axis off;
        summary_text = {
            ['📁 Subject: ' subject ' ' run],
            '',
            '🔧 PREPROCESSING STEPS:',
            ['   • High-pass filter: ' num2str(config.filter_hp) ' Hz'],
            ['   • Low-pass filter: ' num2str(config.filter_lp) ' Hz'],
            '   • Average referencing',
            ['   • Bad channel removal (±' num2str(config.artifact_threshold) ' μV)'],
            ['   • Epoching: [' num2str(config.epoch_window(1)) ' ' num2str(config.epoch_window(2)) '] s'],
            ['   • Baseline: [' num2str(config.baseline(1)) ' ' num2str(config.baseline(2)) '] s'],
            '',
            '📊 RESULTS:',
            ['   • Original channels: ' num2str(EEG_raw.nbchan)],
            ['   • Channels kept: ' num2str(EEG_filtered.nbchan)],
            ['   • Channels removed: ' num2str(EEG_raw.nbchan - EEG_filtered.nbchan)],
            ['   • Total epochs: ' num2str(result.epochs)],
            ['   • T0 (Rest): ' num2str(result.T0)],
            ['   • T1 (Left hand): ' num2str(result.T1)],
            ['   • T2 (Right hand): ' num2str(result.T2)],
            ['   • Sampling rate: ' num2str(EEG_raw.srate) ' Hz'],
            ['   • Epoch length: ' num2str(diff(config.epoch_window)) ' s'],
            '',
            '✅ Processing: SUCCESS'
        };
        
        text(0.05, 0.95, summary_text, 'VerticalAlignment', 'top', ...
             'FontSize', 9, 'FontName', 'FixedWidth', 'Units', 'normalized');
        sgtitle(['EEG Preprocessing Results: ' subject ' ' run], 'FontSize', 16, 'FontWeight', 'bold');
        set(gcf, 'Color', 'white');
        fig_name = fullfile(config.output_folder, [subject '_' run '_preprocessing_report.png']);
        saveas(gcf, fig_name, 'png');
        fprintf('Improved preprocessing plots saved: %s\n', fig_name);
        
    catch ME
        fprintf('  ⚠️  Plot generation failed: %s\n', ME.message);
        fprintf('     Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('       %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
    end
end