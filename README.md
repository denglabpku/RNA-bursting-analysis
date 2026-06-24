# RNA bursting analysis

This folder contains MATLAB scripts used for Sox2 RNA transcriptional-burst detection, tracking, and downstream burst-state modeling for Fig. 5.

The workflow has three steps:

1. Step0_Generate_probability_mask.m: convert Ilastik probability-map H5 files into binary TIFF masks used by Step 1.
2. Step1_RNA_detection.m: detect RNA spots from max-projection images, track one spot per ROI, and quantify RNA intensity.
3. Step2_HMM_modeling.m: merge processed tracks from control and JQ1-treated datasets, infer ON/OFF transcriptional states, and summarize burst statistics.

## Data Organization

The scripts are configured for the Fig. 5 supplementary imaging data folder:

    D:\ImageData\Supplementary_Imaging_Data\Fig5_Sox2RNA_burst_modeling

Expected structure:

    Fig5_Sox2RNA_burst_modeling/
    |-- Untreated/
    |   |-- raw_data/
    |   |-- max_projection/
    |   |-- processed_data/
    |-- JQ1-treated/
    |   |-- raw_data/
    |   |-- max_projection/
    |   |-- processed_data/

File naming is aligned by condition and three-digit index:

- Untreated raw data: RNA-24xMS2_NC_Z_30s_1h_###.tif
- Untreated max projection: RNA-24xMS2_NC_Z_30s_1h_###-Max.tif
- Untreated processed tracks: RNA-24xMS2_NC_Z_30s_1h_###-Max-tracks.mat
- JQ1-treated raw data: RNA-24xMS2_1uMJQ1_Z_30s_1h_###.tif
- JQ1-treated max projection: RNA-24xMS2_1uMJQ1_Z_30s_1h_###-Max.tif
- JQ1-treated processed tracks: RNA-24xMS2_1uMJQ1_Z_30s_1h_###-Max-tracks.mat

Update data_root in Step0, Step1, and Step2 if the data are stored elsewhere.

## Requirements

- MATLAB
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Bioinformatics Toolbox or equivalent HMM functions for hmmviterbi
- Curve Fitting Toolbox
- Helper functions used by the scripts, including TIFFreader, TIFwriter, and slanCM

## Step 0: Generate Probability Masks

Run:

    Step0_Generate_probability_mask

### Inputs

For each condition, Step 0 scans the max_projection folder for files matching:

    *-Max.tif

For every max-projection image, it expects the corresponding Ilastik H5 probability-map file:

    <input-prefix>-Max_Probabilities.h5

The H5 dataset path is:

    /exported_data

### Main Processing Steps

1. Read each max-projection RNA image stack.
2. Read the corresponding Ilastik probability map.
3. Threshold the RNA spot probability channel with probability_threshold = 0.5.
4. Convert the Ilastik export layout back to MATLAB image order.
5. Save the binary probability mask as a TIFF file.
6. Create an empty ROI text file if no ROI file is present.

### Outputs

Step 0 writes intermediate files in the same max_projection folder:

- <input-prefix>-Max-Prob.tif: binary probability mask used by Step 1.
- <input-prefix>-Max-ROI.txt: ROI file placeholder, if it did not already exist.

If the Fig. 5 data package already contains *-Max-Prob.tif and *-Max-ROI.txt files, Step 0 does not need to be rerun unless probability masks are regenerated.

## Step 1: RNA Detection and Tracking

Run:

    Step1_RNA_detection

### Inputs

Step 1 scans both condition-specific max_projection folders for:

    *-Max.tif

For each image, the script expects:

- RNA max-projection time-lapse image stack: *-Max.tif
- Binary probability mask from Step 0: *-Max-Prob.tif
- ROI file defining cell/spot search regions: *-Max-ROI.txt

### Main Processing Steps

1. Load the RNA max-projection image stack and binary probability mask.
2. Identify connected candidate RNA spots in each frame.
3. For each frame, record the brightest pixel in each candidate spot.
4. For each ROI, select the brightest candidate spot within the ROI.
5. Fill missing spot positions by edge carrying and linear interpolation.
6. Fit a bounded 2D Gaussian to a 5-by-5 patch around each tracked spot.
7. Save both Gaussian-fitted intensity and legacy local-background-corrected intensity.
8. Write the resulting track file to the matching processed_data folder.

The Gaussian integral is computed as:

$$
I = 2\pi A\sigma^2
$$

where $A$ is the fitted Gaussian amplitude and $\sigma$ is the fitted Gaussian width.

### Outputs

For each input max-projection image, Step 1 writes:

    <input-prefix>-Max-tracks.mat

The file is saved in the matching processed_data folder.

The saved variables are:

- spots: detected candidate spot coordinates and intensities for each frame.
- tracks: tracked ROI-level RNA spot information.

tracks is a cell array with one row per ROI:

| Field | Content |
| --- | --- |
| tracks{n,1} | Raw detected track: [frame, row, column, max_pixel_intensity] |
| tracks{n,2} | Gap-filled track: [frame, row, column, max_pixel_intensity, gaussian_integral, gaussian_amplitude, gaussian_sigma] |
| tracks{n,3} | Local RNA-window integrated intensity and background |
| tracks{n,4} | Background-corrected RNA-window intensity and background |
| tracks{n,5} | Max-pixel intensity after local background subtraction |

The Gaussian-fitted RNA intensity used downstream is stored in tracks{n,2}(:,5).

## Step 2: HMM Modeling and Burst Analysis

Run:

    Step2_HMM_modeling

### Inputs

Step 2 reads all *-tracks.mat files from the processed_data folders listed in group_dirs.

The current condition labels are:

    group_names = {'Control', 'JQ1 treatment'};
    group_num = [218, 196];
    frame_interval_s = 30;

Update group_num if the number of ROIs per group changes.

### Gaussian Intensity Scaling

Step 2 replaces the legacy max-pixel intensity trace in tracks{n,5} with the Gaussian-fitted intensity when available:

    analysis_intensity(valid_idx) = gaussian_intensity(valid_idx) / gaussian_scale_factor;

The default scaling factor is:

    gaussian_scale_factor = 5.473;

This keeps the Gaussian-fitted intensity on a similar scale to the earlier max-pixel intensity traces used for HMM thresholding.

### Main Processing Steps

1. Merge all track files from untreated and JQ1-treated groups.
2. Build RNA intensity heatmaps across all ROIs and frames.
3. Compare legacy max-pixel intensity with Gaussian-fitted intensity.
4. Inspect intensity distributions and fit a three-component Gaussian mixture model.
5. Infer transcriptional ON/OFF states using a two-state HMM and Viterbi decoding.
6. Remove weak ON calls whose segment maximum is below min_intensity_threshold.
7. Summarize each ON/OFF segment as [state, duration_in_frames, mean_intensity, max_intensity].
8. Compare ON-state intensity distributions between conditions.
9. Quantify burst amplitude, duration, burst count, and burst size.
10. Cluster ROIs by total ON-time length and visualize reordered heatmaps.
11. Extract ON/OFF dwell times and fit one-, two-, and three-exponential models to 1-CDF curves.

### HMM Parameters

The current script uses a fixed two-state model:

    transitionMatrix = [0.8, 0.2;
                        0.2, 0.8];

The emission probabilities are initialized from Gaussian probability density functions over discretized RNA intensity values. HMM states are decoded with hmmviterbi.

### Outputs

Step 2 generates MATLAB figures for:

- Merged RNA-intensity heatmaps
- Max-pixel versus Gaussian-fitted intensity comparison
- RNA intensity histograms and Gaussian mixture fitting
- Example HMM ON/OFF state traces
- ON-state intensity distributions
- Burst amplitude, duration, burst size, and burst count comparisons
- Reordered heatmaps and cluster composition
- ON/OFF dwell-time 1-CDF curves
- Exponential fits to dwell-time distributions

The script also stores inferred states and burst summaries in memory:

| Field | Content |
| --- | --- |
| merged_tracks{n,5} | Analysis intensity trace, using scaled Gaussian intensity where available |
| merged_tracks{n,6} | Inferred transcriptional state trace, where 0 is OFF and 1 is ON |
| merged_tracks{n,7} | Burst segment summary: [state, duration_in_frames, mean_intensity, max_intensity] |
| intensity | Group-level burst amplitude, duration, count, and size summaries |
| on_off_time | Group-level ON and OFF dwell-time distributions |
| param_fitresult | Parameters from exponential fitting of dwell-time 1-CDF curves |

## Notes

- Step 0 is only required when starting from Ilastik H5 probability maps.
- Step 1 should be run before Step 2 if new max-projection image data need to be processed.
- If processed_data already contains *-tracks.mat files, Step 2 can be run directly.
- The scripts are written for two groups: untreated control and JQ1-treated cells.
- The default imaging interval is 30 s per frame.
- The default weak-burst filter in Step 2 is min_intensity_threshold = 1600.

## Citation

If you use this code, please cite the associated paper:

Coordinated dynamics of condensates and enhancer-promoter looping revealed by live-cell imaging.
