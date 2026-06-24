%  Step0_Generate_probability_mask.m
%  Bo Wang, July 2025
close all; clear; clc;

%  DESCRIPTION
%  This script converts Ilastik probability-map outputs into binary TIFF masks
%  used by Step1_RNA_detection.m. It also creates an empty ROI text file when
%  the corresponding ROI file is missing, so users can fill ROI coordinates
%  before running Step 1.

%% Configuration
data_root = 'D:\ImageData\Supplementary_Imaging_Data\Fig5_Sox2RNA_burst_modeling';

group_names = {'Untreated', 'JQ1-treated'};
max_projection_dirs = {fullfile(data_root, 'Untreated', 'max_projection'), ...
                       fullfile(data_root, 'JQ1-treated', 'max_projection')};

probability_threshold = 0.5;

%% Convert Ilastik probability maps to binary masks
for group_iter = 1:length(max_projection_dirs)
    input_dir = max_projection_dirs{group_iter};
    filelist = dir(fullfile(input_dir, '*-Max.tif'));

    fprintf('Processing %s: %d max-projection files\n', group_names{group_iter}, numel(filelist));

    for file_iter = 1:numel(filelist)
        filepath = fullfile(filelist(file_iter).folder, filelist(file_iter).name);
        [folder, name, ~] = fileparts(filepath);
        filepath_no_ext = fullfile(folder, name);

        h5_filepath = [filepath_no_ext, '_Probabilities.h5'];
        prob_tif_path = [filepath_no_ext, '-Prob.tif'];
        roi_filepath = [filepath_no_ext, '-ROI.txt'];

        if ~isfile(h5_filepath)
            warning('Missing Ilastik probability file: %s', h5_filepath);
            continue;
        end

        img_series = TIFFreader(filepath);
        h = size(img_series, 1);
        w = size(img_series, 2);
        numberOfPages = size(img_series, 3);

        % Ilastik exports probability maps as [channel, X, Y, frame]. The RNA
        % spot channel is thresholded and transposed back to MATLAB image order.
        dataset = h5read(h5_filepath, '/exported_data');
        probability_map = reshape(dataset(1, :, :, :), [w, h, numberOfPages]);
        bw = probability_map >= probability_threshold;

        mask = false(h, w, numberOfPages);
        for frame_iter = 1:numberOfPages
            mask(:, :, frame_iter) = bw(:, :, frame_iter)';
        end

        TIFwriter(uint8(mask), prob_tif_path, 'lzw');

        % Step 1 expects an ROI file with four corner points per ROI. Create an
        % empty placeholder only when no ROI file exists yet.
        if ~isfile(roi_filepath)
            fileID = fopen(roi_filepath, 'w');
            if fileID == -1
                warning('Unable to create ROI file: %s', roi_filepath);
            else
                fclose(fileID);
            end
        end

        fprintf('  Wrote %s\n', prob_tif_path);
    end
end
