%  Step2_HMM_modeling.m
%  Bo Wang, July 2025
clc; close all; clear;

%  DESCRIPTION
%  This script merges RNA bursting tracks from different treatment conditions,
%  replaces the legacy max-pixel intensity with Gaussian-fitted RNA intensity
%  when available, and performs HMM-based burst-state calling. The downstream
%  sections summarize burst amplitude, duration, burst size, heatmap clustering,
%  and ON/OFF dwell-time distributions.

%% Configuration
group_dirs = {'D:\ImageData\Supplementary_Imaging_Data\Fig5_Sox2RNA_burst_modeling\Untreated\processed_data\', ...
              'D:\ImageData\Supplementary_Imaging_Data\Fig5_Sox2RNA_burst_modeling\JQ1-treated\processed_data\'};
group_names = {'Control', 'JQ1 treatment'};
group_num = [218, 196];       % Number of ROIs in each group.
frame_interval_s = 30;        % Imaging interval.

% Gaussian integrals are divided by this factor so HMM thresholds remain on
% the same approximate scale as the earlier max-pixel traces.
gaussian_scale_factor = 5.473;

%% Merge track files from all groups
roi_idx = [];
merged_tracks = cell(0);

for group_iter = 1:length(group_dirs)
    dirpath = group_dirs{group_iter};
    filelist = dir([dirpath, '*-tracks.mat']);
    for file_iter = 1:size(filelist, 1)
        filename = filelist(file_iter).name;
        load([dirpath, filename], "tracks");
        merged_tracks = vertcat(merged_tracks, tracks);
        roi_idx = [roi_idx; size(tracks, 1), group_iter];
    end
end

cumsum_num = cumsum(group_num);
numberOfPages = size(merged_tracks{1, 1}, 1);
roi_idx(:, 3) = cumsum(roi_idx(:, 1));

if size(merged_tracks, 1) ~= cumsum_num(end)
    warning('group_num does not match the number of merged tracks. Update group_num before group-level analysis.');
end

group_ranges = cell(length(group_num), 1);
start_idx = 1;
for group_iter = 1:length(group_num)
    group_ranges{group_iter} = start_idx:cumsum_num(group_iter);
    start_idx = cumsum_num(group_iter) + 1;
end

%% Replace max-pixel intensity with Gaussian-fitted RNA intensity
max_intensity_heatmap = zeros(size(merged_tracks, 1), numberOfPages);
gaussian_raw_heatmap = NaN(size(merged_tracks, 1), numberOfPages);
rna_intensity_heatmap = zeros(size(merged_tracks, 1), numberOfPages);

for track_iter = 1:size(merged_tracks, 1)
    l = length(merged_tracks{track_iter, 5});
    max_intensity_heatmap(track_iter, 1:l) = merged_tracks{track_iter, 5}';

    analysis_intensity = merged_tracks{track_iter, 5};

    if size(merged_tracks{track_iter, 2}, 2) >= 5
        gaussian_intensity = merged_tracks{track_iter, 2}(:, 5);
        l_gauss = length(gaussian_intensity);
        gaussian_raw_heatmap(track_iter, 1:l_gauss) = gaussian_intensity';

        valid_idx = ~isnan(gaussian_intensity);
        analysis_intensity(valid_idx) = gaussian_intensity(valid_idx) / gaussian_scale_factor;
    end

    merged_tracks{track_iter, 5} = analysis_intensity;
    l = length(analysis_intensity);
    rna_intensity_heatmap(track_iter, 1:l) = analysis_intensity';
end

%% Plot merged RNA-intensity heatmap
figure;
h1 = heatmap(rna_intensity_heatmap);
colorused = slanCM(97);
colorused = colorused(end:-1:1, :);
colorused = [imresize(colorused(1:128, :), [85, 3], 'bilinear'); colorused(129:256, :)];
h1.Colormap = colorused;
h1.ColorLimits = [0, 6e3];
h1.GridVisible = 'off';
xlabel('Frame');
ylabel('Single Cell');
for i = 1:numberOfPages
    h1.XDisplayLabels{i, 1} = '';
end
for i = 1:size(merged_tracks, 1)
    h1.YDisplayLabels{i, 1} = '';
end

%% Compare legacy max-pixel intensity with Gaussian-fitted intensity
x = max_intensity_heatmap(:);
y = gaussian_raw_heatmap(:) / gaussian_scale_factor;
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

figure;
scatter(x, y, 10, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.1, ...
    'MarkerEdgeAlpha', 0.1);
hold on;
p = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
text(mean(x), max(y), sprintf('y = %.3f x + %.3f', p(1), p(2)), ...
    'FontSize', 12, 'Color', 'r');
xlabel('Max-pixel intensity minus local background');
ylabel('Gaussian intensity / scale factor');
set(gca, 'FontSize', 12, 'LineWidth', 0.9, 'XColor', 'k', 'YColor', 'k');
box on;
xlim([-2.5e3, 1e4]);
ylim([-2.5e3, 1e4]);

figure;
subplot(2, 1, 1);
histogram(x, 'BinWidth', 150);
xlim([-2e3, 8e3]);
xlabel('Max-pixel intensity minus local background', 'Color', 'k');
ylabel('Count', 'Color', 'k');
set(gca, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 0.9);
subplot(2, 1, 2);
histogram(y, 'BinWidth', 300);
xlim([-2e3, 8e3]);
xlabel('Gaussian intensity / scale factor', 'Color', 'k');
ylabel('Count', 'Color', 'k');
set(gca, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 0.9);

%% Inspect Gaussian-intensity distributions with a Gaussian mixture model
figure;
histogram(rna_intensity_heatmap(:), 'BinWidth', 300);
xlabel('Intensity (A.U.)');
ylabel('Count');
xlim([-1500, 8e3]);

figure;
for group_iter = 1:length(group_ranges)
    subplot(length(group_ranges), 1, group_iter);
    histogram(rna_intensity_heatmap(group_ranges{group_iter}, :), 'BinWidth', 300);
    xlabel('Intensity (A.U.)');
    ylabel('Count');
    xlim([-1500, 8e3]);
    ylim([0, 6500]);
    title(group_names{group_iter});
end

gmm_group_idx = 1;
data = rna_intensity_heatmap(group_ranges{gmm_group_idx}, :);
data = data(:);
data = data(~isnan(data));
gm = fitgmdist(data, 3, 'RegularizationValue', 1e-6);

figure;
histogram(data, 'Normalization', 'pdf');
hold on;
x_grid = -1500:10:8e3;
y_total = pdf(gm, x_grid(:));
plot(x_grid, y_total, 'k-', 'LineWidth', 2, 'DisplayName', 'Total GMM');
for i = 1:gm.NumComponents
    mu = gm.mu(i);
    sigma = sqrt(gm.Sigma(:, :, i));
    weight = gm.ComponentProportion(i);
    y_i = weight * normpdf(x_grid, mu, sigma);
    plot(x_grid, y_i, '--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Component %d', i));
end
xlim([-1500, 8e3]);
title([group_names{gmm_group_idx}, ' Gaussian Mixture Model']);
xlabel('Intensity (A.U.)');
ylabel('Probability Density');
legend show;

%% Hidden Markov model
min_intensity = 0;
max_intensity = 3000;
min_intensity_threshold = 1600;

for track_iter = 1:size(merged_tracks, 1)
    observ = rna_intensity_heatmap(track_iter, :);
    observ(observ < min_intensity) = min_intensity;
    observ(observ > max_intensity) = max_intensity;
    observ = round(observ) + 1;

    transitionMatrix = [0.8, 0.2;
                        0.2, 0.8];
    emissionMatrix = [normpdf(0:1:max_intensity, max_intensity/8, max_intensity/8);
                      normpdf(0:1:max_intensity, max_intensity/3*2, max_intensity/3)];

    estimatedStates = hmmviterbi(observ, transitionMatrix, emissionMatrix);
    merged_tracks{track_iter, 6} = (estimatedStates - 1)';

    % Summarize consecutive OFF/ON segments as
    % [state, duration_in_frames, mean_intensity, max_intensity].
    rna_state = merged_tracks{track_iter, 6};
    changes = [0, find(diff(rna_state))', length(rna_state)];
    lengths = diff(changes);
    values = rna_state(changes(2:end));
    merged_tracks{track_iter, 7} = [values, lengths'];

    start_num = 0;
    for i = 1:size(merged_tracks{track_iter, 7}, 1)
        frame_idx = (start_num+1):(start_num+merged_tracks{track_iter, 7}(i, 2));
        merged_tracks{track_iter, 7}(i, 3) = mean(rna_intensity_heatmap(track_iter, frame_idx));
        merged_tracks{track_iter, 7}(i, 4) = max(rna_intensity_heatmap(track_iter, frame_idx));
        start_num = start_num + merged_tracks{track_iter, 7}(i, 2);
    end

    % Remove weak ON calls whose segment maximum is below the intensity
    % threshold, then rebuild the state vector and segment summary.
    merged_tracks{track_iter, 7}(merged_tracks{track_iter, 7}(:, 4) < min_intensity_threshold, 1) = 0;
    new_state = [];
    for i = 1:size(merged_tracks{track_iter, 7}, 1)
        new_state = [new_state; repmat(merged_tracks{track_iter, 7}(i, 1), ...
            [merged_tracks{track_iter, 7}(i, 2), 1])];
    end
    merged_tracks{track_iter, 6} = new_state;

    rna_state = merged_tracks{track_iter, 6};
    changes = [0, find(diff(rna_state))', length(rna_state)];
    lengths = diff(changes);
    values = rna_state(changes(2:end));
    merged_tracks{track_iter, 7} = [values, lengths'];

    start_num = 0;
    for i = 1:size(merged_tracks{track_iter, 7}, 1)
        frame_idx = (start_num+1):(start_num+merged_tracks{track_iter, 7}(i, 2));
        merged_tracks{track_iter, 7}(i, 3) = mean(rna_intensity_heatmap(track_iter, frame_idx));
        merged_tracks{track_iter, 7}(i, 4) = max(rna_intensity_heatmap(track_iter, frame_idx));
        start_num = start_num + merged_tracks{track_iter, 7}(i, 2);
    end
end

%% Compare active-state intensity at one frame and across all ON frames
frame_idx = 50;
if frame_idx <= numberOfPages && length(group_ranges) >= 2
    frame_on_intensity = cell(length(group_ranges), 1);
    for group_iter = 1:length(group_ranges)
        for track_iter = group_ranges{group_iter}
            if merged_tracks{track_iter, 6}(frame_idx) == 1
                frame_on_intensity{group_iter} = [frame_on_intensity{group_iter}; ...
                    merged_tracks{track_iter, 5}(frame_idx)];
            end
        end
        frame_on_intensity{group_iter} = frame_on_intensity{group_iter}(frame_on_intensity{group_iter} > 600);
        fprintf('%s mean active intensity at frame %d: %.3f\n', ...
            group_names{group_iter}, frame_idx, mean(frame_on_intensity{group_iter}));
    end
end

on_state_intensity = cell(length(group_ranges), 1);
for group_iter = 1:length(group_ranges)
    for track_iter = group_ranges{group_iter}
        on_state_intensity{group_iter} = [on_state_intensity{group_iter}; ...
            merged_tracks{track_iter, 5}(merged_tracks{track_iter, 6} == 1)];
    end
end

figure;
subplot(1, 2, 1);
hold on;
for group_iter = 1:length(group_ranges)
    histogram(on_state_intensity{group_iter}(on_state_intensity{group_iter} >= 600), ...
        'BinWidth', 300, 'Normalization', 'probability');
end
hold off;
xlabel('Intensity (A.U.)');
ylabel('Probability');
xlim([0, 9e3]);
legend(group_names);

subplot(1, 2, 2);
hold on;
for group_iter = 1:length(group_ranges)
    histogram(on_state_intensity{group_iter}(on_state_intensity{group_iter} >= 600), ...
        'BinWidth', 300);
end
hold off;
xlabel('Intensity (A.U.)');
ylabel('Counts');
xlim([0, 9e3]);
legend(group_names);

%% Check HMM calls for specified and random ROIs
example_idx = [80, min(326, size(merged_tracks, 1))];
example_num = length(example_idx);
figure;
for j = 1:example_num
    track_iter = example_idx(j);
    vcount = merged_tracks{track_iter, 7};
    vcount(:, 3) = vcount(:, 1) .* vcount(:, 3);
    output = [];
    for i = 1:size(vcount, 1)
        output = [output; repmat(vcount(i, 3), vcount(i, 2), 1)];
    end
    subplot(example_num, 1, j);
    plot(rna_intensity_heatmap(track_iter, :) * gaussian_scale_factor);
    hold on;
    plot(output * gaussian_scale_factor);
    hold off;
    xlim([1, min(120, numberOfPages)]);
    ylim([-1e4, 4e4]);
end

example_num = min(36, length(merged_tracks));
rng(123);
example_idx = randsample(length(merged_tracks), example_num);
example_cols = 3;
example_rows = ceil(example_num / example_cols);
figure;
for j = 1:example_num
    track_iter = example_idx(j);
    vcount = merged_tracks{track_iter, 7};
    vcount(:, 3) = vcount(:, 1) .* vcount(:, 3);
    output = [];
    for i = 1:size(vcount, 1)
        output = [output; repmat(vcount(i, 3), vcount(i, 2), 1)];
    end
    subplot(example_rows, example_cols, j);
    plot(rna_intensity_heatmap(track_iter, :));
    hold on;
    plot(output);
    hold off;
    xlim([1, numberOfPages]);
    ylim([-1e3, 4e3]);
end

%% Burst amplitude, duration, counts, and size comparison
intensity = cell(length(group_ranges), 6);
for group_iter = 1:length(group_ranges)
    temp_intensity = [];
    temp_length = [];
    temp_num = [];
    temp_total_length = [];
    temp_total_size = [];

    for track_iter = group_ranges{group_iter}
        vcount = merged_tracks{track_iter, 7};
        vcount(:, 3) = vcount(:, 1) .* vcount(:, 3);
        temp_intensity = [temp_intensity; vcount(vcount(:, 1) == 1, 3)];
        temp_length = [temp_length; vcount(vcount(:, 1) == 1, 2)];
        temp_num = [temp_num; sum(vcount(:, 1) == 1)];
        temp_total_length = [temp_total_length; sum(vcount(vcount(:, 1) == 1, 2))];
        temp_total_size = [temp_total_size; sum(vcount(vcount(:, 1) == 1, 3) .* vcount(vcount(:, 1) == 1, 2))];
    end

    valid_burst = temp_intensity < 1e4;
    intensity{group_iter, 1} = temp_intensity(valid_burst) * gaussian_scale_factor; % mean amplitude
    intensity{group_iter, 2} = temp_length(valid_burst);                            % duration in frames
    intensity{group_iter, 3} = temp_num;                                             % burst counts
    intensity{group_iter, 4} = intensity{group_iter, 1} .* intensity{group_iter, 2};  % burst size
    intensity{group_iter, 5} = temp_total_length;                                    % total ON duration
    intensity{group_iter, 6} = temp_total_size * gaussian_scale_factor;              % total burst size
end

figure;
subplot(1, 4, 1);
hold on;
for group_iter = 1:length(group_ranges)
    histogram(intensity{group_iter, 1}, 'BinWidth', 2000);
end
hold off;
xlim([0, 40000]);
xlabel('Mean amplitude');

subplot(1, 4, 2);
hold on;
for group_iter = 1:length(group_ranges)
    histogram(intensity{group_iter, 2} * frame_interval_s / 60, 'BinWidth', 2);
end
hold off;
xlim([-5, 70]);
xlabel('Duration (min)');

subplot(1, 4, 3);
hold on;
for group_iter = 1:length(group_ranges)
    histogram(intensity{group_iter, 4}, 'BinWidth', 1.8e4);
end
hold off;
xlim([-2e4, 3e5]);
xlabel('Burst size');

subplot(1, 4, 4);
hold on;
for group_iter = 1:length(group_ranges)
    histogram(intensity{group_iter, 3});
end
hold off;
xlim([-1, 12]);
xlabel('Burst count per ROI');
legend(group_names);

color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];
alpha = 0.4;
figure;
scatter(intensity{1, 1}, intensity{1, 2} * frame_interval_s / 60, 30, ...
    'MarkerFaceColor', color1, 'MarkerEdgeColor', color1, ...
    'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha, 'DisplayName', group_names{1});
hold on;
scatter(intensity{2, 1}, intensity{2, 2} * frame_interval_s / 60, 30, ...
    'MarkerFaceColor', color2, 'MarkerEdgeColor', color2, ...
    'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha, 'DisplayName', group_names{2});
hold off;
xlabel('Intensity (A.U.)');
ylabel('Duration (min)');
legend show;
title('Burst intensity vs duration');
xlim([0, 4e4]);
ylim([0, 65]);

%% Cluster ROIs by total ON-time length and visualize reordered heatmaps
on_time_length_intensity = zeros(size(merged_tracks, 1), 4);
for i = 1:size(merged_tracks, 1)
    on_time = sum(merged_tracks{i, 6});
    total_intensity = sum(merged_tracks{i, 6}(1:length(merged_tracks{i, 5})) .* merged_tracks{i, 5});
    on_time_length_intensity(i, 1:3) = [i, on_time, total_intensity];
    if on_time > 0
        on_time_length_intensity(i, 4) = total_intensity / on_time;
    else
        on_time_length_intensity(i, 4) = NaN;
    end
end

figure;
histogram(on_time_length_intensity(:, 2), 'BinWidth', 3);
xlabel('ON-time length (frames)');
ylabel('Count');

sep_point = [0, 25, 46, 73, numberOfPages];
max_intensity_filter = 6e3;
rna_intensity_reorder = rna_intensity_heatmap;
rna_intensity_reorder(rna_intensity_reorder > max_intensity_filter) = max_intensity_filter;
rna_intensity_reorder(rna_intensity_reorder < 0) = 0;

total_idx = [];
for cluster_iter = 1:length(sep_point)-1
    idx = on_time_length_intensity( ...
        on_time_length_intensity(:, 2) > sep_point(cluster_iter) & ...
        on_time_length_intensity(:, 2) <= sep_point(cluster_iter+1), 1);

    if isempty(idx)
        continue;
    elseif isscalar(idx)
        total_idx = [total_idx; [idx, cluster_iter]];
        continue;
    end

    temp_intensity = rna_intensity_reorder(idx, :);
    cgo = clustergram(temp_intensity, 'Cluster', 'column');
    set(cgo, 'Linkage', 'complete');
    row_order = str2double(cgo.RowLabels);
    idx = idx(row_order);
    idx = idx(end:-1:1);
    total_idx = [total_idx; [idx, ones(size(idx))*cluster_iter]];
end
total_idx = total_idx(end:-1:1, :);

rna_heatmap_reorder = rna_intensity_heatmap(total_idx(:, 1), :);
figure;
h = heatmap(rna_heatmap_reorder);
colorused = slanCM(97);
colorused = colorused(end:-1:1, :);
colorused = [imresize(colorused(1:128, :), [85, 3], 'bilinear'); colorused(129:256, :)];
h.Colormap = colorused;
h.ColorLimits = [0, 6e3];
h.GridVisible = 'off';
xlabel('Frame');
for i = 1:numberOfPages
    h.XDisplayLabels{i, 1} = '';
end

total_idx(:, 3) = 1;
for i = 1:length(group_num)-1
    total_idx(total_idx(:, 1) > cumsum_num(i), 3) = i + 1;
end

cluster_num = zeros(length(sep_point)-1, length(group_num));
for i = 1:length(sep_point)-1
    for j = 1:length(group_num)
        cluster_num(i, j) = sum(total_idx(:, 2) == i & total_idx(:, 3) == j) / ...
            sum(total_idx(:, 3) == j);
    end
end
cluster_labels = arrayfun(@(x) sprintf('Cluster %d', x), ...
    1:(length(sep_point)-1), 'UniformOutput', false);
figure;
heatmap(group_names, cluster_labels, cluster_num);

%% Plot one group heatmap with cluster labels
selected_group_idx = 1; % 1 for Untreated; 2 for JQ1-treated
selected_rows = total_idx(total_idx(:, 3) == selected_group_idx, :);

rna_heatmap_reorder = rna_intensity_heatmap(selected_rows(:, 1), :) * gaussian_scale_factor;
figure;
h = heatmap(rna_heatmap_reorder);
colorused = slanCM(97);
colorused = colorused(end:-1:1, :);
colorused = [imresize(colorused(1:128, :), [85, 3], 'bilinear'); colorused(129:256, :)];
h.Colormap = colorused;
h.ColorLimits = [0, 3e4];
h.GridVisible = 'off';
xlabel('Frame');
title([group_names{selected_group_idx}, ' reordered Gaussian-intensity heatmap']);
for i = 1:numberOfPages
    h.XDisplayLabels{i, 1} = '';
end
for i = 1:size(selected_rows, 1)
    h.YDisplayLabels{i, 1} = '';
end

%% ON/OFF dwell-time distributions
on_off_time = cell(length(group_ranges), 2);
for group_iter = 1:length(group_ranges)
    temp_on_off_time = [];
    for track_iter = group_ranges{group_iter}
        temp_on_off_time = [temp_on_off_time; merged_tracks{track_iter, 7}];
    end

    on_time = temp_on_off_time(temp_on_off_time(:, 1) == 1, 2) * frame_interval_s;
    off_time = temp_on_off_time(temp_on_off_time(:, 1) == 0, 2) * frame_interval_s;
    on_off_time{group_iter, 1} = on_time;
    on_off_time{group_iter, 2} = off_time;

    fprintf('%s mean ON time: %.3f s\n', group_names{group_iter}, mean(on_time));
    fprintf('%s mean OFF time: %.3f s\n', group_names{group_iter}, mean(off_time));
end

group_idx = 1;
c1 = cdfplot(on_off_time{group_idx, 1}/60);
on_x = c1.XData;
on_y = 1 - c1.YData;
close;
c2 = cdfplot(on_off_time{group_idx, 2}/60);
off_x = c2.XData;
off_y = 1 - c2.YData;
close;
figure;
plot(on_x, on_y);
hold on;
plot(off_x, off_y);
hold off;
legend({'ON time', 'OFF time'});
ylabel('1-CDF');
xlabel('Time (min)');
xlim([0, 30]);
title(group_names{group_idx});

%% Fit one-, two-, and three-exponential models to ON/OFF 1-CDF curves
group_idx = 1;  % 1: control; 2: JQ1 treatment
on_off_idx = 1; % 1: ON time; 2: OFF time
on_off_labels = {'ON Time (s)', 'OFF Time (s)'};
[y, x] = ecdf(on_off_time{group_idx, on_off_idx});
y = 1 - y;
x(1) = 0;
param_fitresult = zeros(9, 3);

figure;
exp1 = fittype('exp(a*x)');
[f, gof, output] = fit(x, y, exp1, ...
    'StartPoint', -0.01, ...
    'Lower', -1, ...
    'Upper', 0);
subplot(2, 3, 1);
plot(f, x, y);
xlabel(on_off_labels{on_off_idx});
ylabel('1-CDF');
title('1-exp');
param_fitresult(1:3, 1) = [1, -f.a, -1/f.a];

subplot(2, 3, 4);
plot(x, output.residuals);
xlabel(on_off_labels{on_off_idx});
ylabel('Residuals');
title(['R-square: ', num2str(gof.rsquare)]);
ylim([-0.05, 0.05]);

exp2 = fittype('a*exp(b*x) + (1-a)*exp(d*x)');
[f, gof, output] = fit(x, y, exp2, ...
    'StartPoint', [0.5, -0.005, -0.01], ...
    'Lower', [0, -1, -1], ...
    'Upper', [1, 0, 0]);
subplot(2, 3, 2);
plot(f, x, y);
xlabel(on_off_labels{on_off_idx});
ylabel('1-CDF');
title('2-exp');
param_fitresult(1:6, 2) = [f.a, -f.b, -1/f.b, 1-f.a, -f.d, -1/f.d];

subplot(2, 3, 5);
plot(x, output.residuals);
xlabel(on_off_labels{on_off_idx});
ylabel('Residuals');
title(['R-square: ', num2str(gof.rsquare)]);
ylim([-0.05, 0.05]);

exp3 = fittype('a*exp(b*x) + c*exp(d*x) + (1-a-c)*exp(f*x)');
[f, gof, output] = fit(x, y, exp3, ...
    'StartPoint', [1, -0.005, 1, -0.01, -0.05], ...
    'Lower', [0, -1, 0, -1, -1], ...
    'Upper', [1, 0, 1, 0, 0]);
subplot(2, 3, 3);
plot(f, x, y);
xlabel(on_off_labels{on_off_idx});
ylabel('1-CDF');
title('3-exp');
param_fitresult(1:9, 3) = [f.a, -f.b, -1/f.b, f.c, -f.d, -1/f.d, ...
    1-f.a-f.c, -f.f, -1/f.f];

subplot(2, 3, 6);
plot(x, output.residuals);
xlabel(on_off_labels{on_off_idx});
ylabel('Residuals');
title(['R-square: ', num2str(gof.rsquare)]);
ylim([-0.05, 0.05]);

fprintf('%s mean %s: %.3f s\n', group_names{group_idx}, ...
    on_off_labels{on_off_idx}, mean(on_off_time{group_idx, on_off_idx}));
