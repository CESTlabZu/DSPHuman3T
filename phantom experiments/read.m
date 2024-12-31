% Clear the workspace, close all figure windows, and clear the command window
clear all;
close all;
clc;

%% Define matrix sizes and parameters
nx = 144; % Number of columns in the matrix
ny = 144; % Number of rows in the matrix

wl = 0.5; % Low saturation power
wh = 1;   % High saturation power

% Initialize matrices for storing results
Image_MTR_2pool = zeros(nx, ny, 47, 8);
Image_AREX_2pool = zeros(nx, ny, 47, 8);

% RF frequency offset in Hz
x = [-10 -8 -6 -5 -4.75 -4.5 -4.25 -4 -3.75 -3.5 -3.25 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 6 8 10]' * 127.6596;

% RF frequency offset for LD fitting in Hz
x_2pool = [-10 -8 -0.5 -0.25 0 0.25 0.5 8 10]' * 127.6596;

% RF frequency array and other parameters
FreqArray = [-10, -8, -6, -5:0.25:5, 6, 8, 10]; % RF frequency array
deta_cs = 0.025;                               % Step size for RF frequency
range_cs = -10:deta_cs:10;                     % Range for RF frequency correction
interp1NOE = 801;                              % Interpolation parameter
referenceCriticalValue = 100;                  % Reference critical value
cest1_0 = 0;                                   % CEST parameter

%% Load data
load('data\roi_brain');

%% Read and normalize CEST data
for ii = 1:51
    nImag1_0p5uT(:, :, ii) = double(Imag1_0p5uT(:, :, ii)) ./ mean(Imag1_0p5uT(:, :, [1, 2, 50, 51]), 3);
    nImag1_1uT(:, :, ii) = double(Imag1_1p0uT(:, :, ii)) ./ mean(Imag1_1p0uT(:, :, [1, 2, 50, 51]), 3);
end

%% B0 correction for normalized images
[nImag1_0p5uT_corr] = B0Correction(nImag1_0p5uT, roi_whole, FreqArray, range_cs, interp1NOE, deta_cs);
[nImag1_1p0uT_corr] = B0Correction(nImag1_1uT, roi_whole, FreqArray, range_cs, interp1NOE, deta_cs);

%% Calculate DSP images
nImag1_DSP_0p5uT_corr = 1 ./ (1 + (1 ./ nImag1_1p0uT_corr - 1) .* ((wl / wh)^2));
Image_MTR_APT_0p5uT_corr = -(squeeze(nImag1_0p5uT_corr(:, :,38) - nImag1_DSP_0p5uT_corr(:, :, 38)) - (squeeze(nImag1_0p5uT_corr(:, :, 44) - nImag1_DSP_0p5uT_corr(:, :, 44)))) .* roi_whole;
Image_MTR_NOE_0p5uT_corr = -(squeeze(nImag1_0p5uT_corr(:, :, 10) - nImag1_DSP_0p5uT_corr(:, :, 10)) - (squeeze(nImag1_0p5uT_corr(:, :, 44) - nImag1_DSP_0p5uT_corr(:, :, 44)))) .* roi_whole;

Image_AREX_APT_0p5uT_corr = (squeeze(1 ./ nImag1_0p5uT_corr(:, :, 38) - 1 ./ nImag1_DSP_0p5uT_corr(:, :, 38)) - (squeeze(1 ./ nImag1_0p5uT_corr(:, :, 44) - 1 ./ nImag1_DSP_0p5uT_corr(:, :, 44)))) .* 1000 ./ double(Imag1_T1) .* roi_whole;
Image_AREX_NOE_0p5uT_corr = (squeeze(1 ./ nImag1_0p5uT_corr(:, :, 10) - 1 ./ nImag1_DSP_0p5uT_corr(:, :, 10)) - (squeeze(1 ./ nImag1_0p5uT_corr(:, :, 44) - 1 ./ nImag1_DSP_0p5uT_corr(:, :, 44)))) .* 1000 ./ double(Imag1_T1) .* roi_whole;

%% Calculate Z-spectra for different phantoms
for ii = 1:47
    Zspectra_0p5uT_eggwhite_corr(ii) = nansum(nansum(nImag1_0p5uT_corr(:, :, ii) .* roi_eggwhite)) / nansum(nansum(roi_eggwhite));
    Zspectra_0p5uT_glu_pH7p2_corr(ii) = nansum(nansum(nImag1_0p5uT_corr(:, :, ii) .* roi_glu_pH7p2)) / nansum(nansum(roi_glu_pH7p2));
    Zspectra_0p5uT_glu_pH7p0_corr(ii) = nansum(nansum(nImag1_0p5uT_corr(:, :, ii) .* roi_glu_pH7p0)) / nansum(nansum(roi_glu_pH7p0));
    Zspectra_0p5uT_glu_pH6p5_corr(ii) = nansum(nansum(nImag1_0p5uT_corr(:, :, ii) .* roi_glu_pH6p5)) / nansum(nansum(roi_glu_pH6p5));
    Zspectra_0p5uT_MnCl2_0p04mM_corr(ii) = nansum(nansum(nImag1_0p5uT_corr(:, :, ii) .* roi_MnCl2_0p04mM)) / nansum(nansum(roi_MnCl2_0p04mM));
    Zspectra_0p5uT_MnCl2_0p08mM_corr(ii) = nansum(nansum(nImag1_0p5uT_corr(:, :, ii) .* roi_MnCl2_0p08mM)) / nansum(nansum(roi_MnCl2_0p08mM));

    Zspectra_1p0uT_eggwhite_corr(ii) = nansum(nansum(nImag1_1p0uT_corr(:, :, ii) .* roi_eggwhite)) / nansum(nansum(roi_eggwhite));
    Zspectra_1p0uT_glu_pH7p2_corr(ii) = nansum(nansum(nImag1_1p0uT_corr(:, :, ii) .* roi_glu_pH7p2)) / nansum(nansum(roi_glu_pH7p2));
    Zspectra_1p0uT_glu_pH7p0_corr(ii) = nansum(nansum(nImag1_1p0uT_corr(:, :, ii) .* roi_glu_pH7p0)) / nansum(nansum(roi_glu_pH7p0));
    Zspectra_1p0uT_glu_pH6p5_corr(ii) = nansum(nansum(nImag1_1p0uT_corr(:, :, ii) .* roi_glu_pH6p5)) / nansum(nansum(roi_glu_pH6p5));
    Zspectra_1p0uT_MnCl2_0p04mM_corr(ii) = nansum(nansum(nImag1_1p0uT_corr(:, :, ii) .* roi_MnCl2_0p04mM)) / nansum(nansum(roi_MnCl2_0p04mM));
    Zspectra_1p0uT_MnCl2_0p08mM_corr(ii) = nansum(nansum(nImag1_1p0uT_corr(:, :, ii) .* roi_MnCl2_0p08mM)) / nansum(nansum(roi_MnCl2_0p08mM));
end

%% Smooth the Z-spectra data
 d_Zspectra_0p5uT_eggwhite_corr=Zspectra_0p5uT_eggwhite_corr;
 d_Zspectra_0p5uT_glu_pH7p2_corr=Zspectra_0p5uT_glu_pH7p2_corr;
 d_Zspectra_0p5uT_glu_pH7p0_corr=Zspectra_0p5uT_glu_pH7p0_corr;
 d_Zspectra_0p5uT_glu_pH6p5_corr=Zspectra_0p5uT_glu_pH6p5_corr;
 d_Zspectra_0p5uT_MnCl2_0p04mM_corr=Zspectra_0p5uT_MnCl2_0p04mM_corr;
 d_Zspectra_0p5uT_MnCl2_0p08mM_corr=Zspectra_0p5uT_MnCl2_0p08mM_corr;

 d_Zspectra_1p0uT_eggwhite_corr=Zspectra_1p0uT_eggwhite_corr;
 d_Zspectra_1p0uT_glu_pH7p2_corr=Zspectra_1p0uT_glu_pH7p2_corr;
 d_Zspectra_1p0uT_glu_pH7p0_corr=Zspectra_1p0uT_glu_pH7p0_corr;
 d_Zspectra_1p0uT_glu_pH6p5_corr=Zspectra_1p0uT_glu_pH6p5_corr;
 d_Zspectra_1p0uT_MnCl2_0p04mM_corr=Zspectra_1p0uT_MnCl2_0p04mM_corr;
 d_Zspectra_1p0uT_MnCl2_0p08mM_corr=Zspectra_1p0uT_MnCl2_0p08mM_corr; 
 
threshold = 4;
smooth_range = 4:44;
d_Zspectra_0p5uT_eggwhite_corr(smooth_range) = smoothdata(Zspectra_0p5uT_eggwhite_corr(smooth_range), "gaussian", threshold);
d_Zspectra_0p5uT_glu_pH7p2_corr(smooth_range) = smoothdata(Zspectra_0p5uT_glu_pH7p2_corr(smooth_range), "gaussian", threshold);
d_Zspectra_0p5uT_glu_pH7p0_corr(smooth_range) = smoothdata(Zspectra_0p5uT_glu_pH7p0_corr(smooth_range), "gaussian", threshold);
d_Zspectra_0p5uT_glu_pH6p5_corr(smooth_range) = smoothdata(Zspectra_0p5uT_glu_pH6p5_corr(smooth_range), "gaussian", threshold);
d_Zspectra_0p5uT_MnCl2_0p04mM_corr(smooth_range) = smoothdata(Zspectra_0p5uT_MnCl2_0p04mM_corr(smooth_range), "gaussian", threshold);
d_Zspectra_0p5uT_MnCl2_0p08mM_corr(smooth_range) = smoothdata(Zspectra_0p5uT_MnCl2_0p08mM_corr(smooth_range), "gaussian", threshold);

d_Zspectra_1p0uT_eggwhite_corr(smooth_range) = smoothdata(Zspectra_1p0uT_eggwhite_corr(smooth_range), "gaussian", threshold);
d_Zspectra_1p0uT_glu_pH7p2_corr(smooth_range) = smoothdata(Zspectra_1p0uT_glu_pH7p2_corr(smooth_range), "gaussian", threshold);
d_Zspectra_1p0uT_glu_pH7p0_corr(smooth_range) = smoothdata(Zspectra_1p0uT_glu_pH7p0_corr(smooth_range), "gaussian", threshold);
d_Zspectra_1p0uT_glu_pH6p5_corr(smooth_range) = smoothdata(Zspectra_1p0uT_glu_pH6p5_corr(smooth_range), "gaussian", threshold);
d_Zspectra_1p0uT_MnCl2_0p04mM_corr(smooth_range) = smoothdata(Zspectra_1p0uT_MnCl2_0p04mM_corr(smooth_range), "gaussian", threshold);
d_Zspectra_1p0uT_MnCl2_0p08mM_corr(smooth_range) = smoothdata(Zspectra_1p0uT_MnCl2_0p08mM_corr(smooth_range), "gaussian", threshold);

%% Calculate DSP Z-spectra for different phantoms
Zspectra_DSP_eggwhite_0p5uT_corr = 1 ./ (1 + (1 ./ d_Zspectra_1p0uT_eggwhite_corr - 1) .* ((wl / wh)^2));
Zspectra_DSP_glu_pH7p2_0p5uT_corr = 1 ./ (1 + (1 ./ d_Zspectra_1p0uT_glu_pH7p2_corr - 1) .* ((wl / wh)^2));
Zspectra_DSP_glu_pH7p0_0p5uT_corr = 1 ./ (1 + (1 ./ d_Zspectra_1p0uT_glu_pH7p0_corr - 1) .* ((wl / wh)^2));
Zspectra_DSP_glu_pH6p5_0p5uT_corr = 1 ./ (1 + (1 ./ d_Zspectra_1p0uT_glu_pH6p5_corr - 1) .* ((wl / wh)^2));
Zspectra_DSP_MnCl2_0p04mM_0p5uT_corr = 1 ./ (1 + (1 ./ d_Zspectra_1p0uT_MnCl2_0p04mM_corr - 1) .* ((wl / wh)^2));
Zspectra_DSP_MnCl2_0p08mM_0p5uT_corr = 1 ./ (1 + (1 ./ d_Zspectra_1p0uT_MnCl2_0p08mM_corr - 1) .* ((wl / wh)^2));

%% Calculate MTR spectra
Zspectra_MTR_eggwhiteT_corr = Zspectra_DSP_eggwhite_0p5uT_corr - d_Zspectra_0p5uT_eggwhite_corr;
Zspectra_MTR_glu_pH7p2T_corr = Zspectra_DSP_glu_pH7p2_0p5uT_corr - d_Zspectra_0p5uT_glu_pH7p2_corr;
Zspectra_MTR_glu_pH7p0T_corr = Zspectra_DSP_glu_pH7p0_0p5uT_corr - d_Zspectra_0p5uT_glu_pH7p0_corr;
Zspectra_MTR_glu_pH6p5T_corr = Zspectra_DSP_glu_pH6p5_0p5uT_corr - d_Zspectra_0p5uT_glu_pH6p5_corr;
Zspectra_MTR_MnCl2_0p04mMT_corr = Zspectra_DSP_MnCl2_0p04mM_0p5uT_corr - d_Zspectra_0p5uT_MnCl2_0p04mM_corr;
Zspectra_MTR_MnCl2_0p08mMT_corr = Zspectra_DSP_MnCl2_0p08mM_0p5uT_corr - d_Zspectra_0p5uT_MnCl2_0p08mM_corr;

% R1 map calculation
R1_map = 1000 ./ double(Imag1_T1);
R1_map(R1_map == inf | R1_map == -inf) = nan;

% R1 values for different phantoms
Value_R1_eggwhite = nansum(nansum(R1_map .* roi_eggwhite)) / nansum(nansum(roi_eggwhite));
Value_R1_glu_pH7p2 = nansum(nansum(R1_map .* roi_glu_pH7p2)) / nansum(nansum(roi_glu_pH7p2));
Value_R1_glu_pH7p0 = nansum(nansum(R1_map .* roi_glu_pH7p0)) / nansum(nansum(roi_glu_pH7p0));
Value_R1_glu_pH6p5 = nansum(nansum(R1_map .* roi_glu_pH6p5)) / nansum(nansum(roi_glu_pH6p5));
Value_R1_MnCl2_0p04mM = nansum(nansum(R1_map .* roi_MnCl2_0p04mM)) / nansum(nansum(roi_MnCl2_0p04mM));
Value_R1_MnCl2_0p08mM = nansum(nansum(R1_map .* roi_MnCl2_0p08mM)) / nansum(nansum(roi_MnCl2_0p08mM));

% Calculate AREX spectra
Zspectra_AREX_eggwhiteT_corr = (1 ./ d_Zspectra_0p5uT_eggwhite_corr - 1 ./ Zspectra_DSP_eggwhite_0p5uT_corr) .* Value_R1_eggwhite;
Zspectra_AREX_glu_pH7p2T_corr = (1 ./ d_Zspectra_0p5uT_glu_pH7p2_corr - 1 ./ Zspectra_DSP_glu_pH7p2_0p5uT_corr) .* Value_R1_glu_pH7p2;
Zspectra_AREX_glu_pH7p0T_corr = (1 ./ d_Zspectra_0p5uT_glu_pH7p0_corr - 1 ./ Zspectra_DSP_glu_pH7p0_0p5uT_corr) .* Value_R1_glu_pH7p0;
Zspectra_AREX_glu_pH6p5T_corr = (1 ./ d_Zspectra_0p5uT_glu_pH6p5_corr - 1 ./ Zspectra_DSP_glu_pH6p5_0p5uT_corr) .* Value_R1_glu_pH6p5;
Zspectra_AREX_MnCl2_0p04mMT_corr = (1 ./ d_Zspectra_0p5uT_MnCl2_0p04mM_corr - 1 ./ Zspectra_DSP_MnCl2_0p04mM_0p5uT_corr) .* Value_R1_MnCl2_0p04mM;
Zspectra_AREX_MnCl2_0p08mMT_corr = (1 ./ d_Zspectra_0p5uT_MnCl2_0p08mM_corr - 1 ./ Zspectra_DSP_MnCl2_0p08mM_0p5uT_corr) .* Value_R1_MnCl2_0p08mM;

%% Calculate MTR asymmetry for different phantoms
for ii = 1:23
    MTRasym_0p5uT_eggwhite_corr(ii) = d_Zspectra_0p5uT_eggwhite_corr(24-ii) - d_Zspectra_0p5uT_eggwhite_corr(24+ii);
    MTRasym_0p5uT_glu_pH7p2_corr(ii) = d_Zspectra_0p5uT_glu_pH7p2_corr(24-ii) - d_Zspectra_0p5uT_glu_pH7p2_corr(24+ii);
    MTRasym_0p5uT_glu_pH7p0_corr(ii) = d_Zspectra_0p5uT_glu_pH7p0_corr(24-ii) - d_Zspectra_0p5uT_glu_pH7p0_corr(24+ii);
    MTRasym_0p5uT_glu_pH6p5_corr(ii) = d_Zspectra_0p5uT_glu_pH6p5_corr(24-ii) - d_Zspectra_0p5uT_glu_pH6p5_corr(24+ii);
    MTRasym_0p5uT_MnCl2_0p04mM_corr(ii) = d_Zspectra_0p5uT_MnCl2_0p04mM_corr(24-ii) - d_Zspectra_0p5uT_MnCl2_0p04mM_corr(24+ii);
    MTRasym_0p5uT_MnCl2_0p08mM_corr(ii) = d_Zspectra_0p5uT_MnCl2_0p08mM_corr(24-ii) - d_Zspectra_0p5uT_MnCl2_0p08mM_corr(24+ii);
end

% Flip the MTR asymmetry data
fMTRasym_0p5uT_eggwhite_corr = flip(MTRasym_0p5uT_eggwhite_corr);
fMTRasym_0p5uT_glu_pH7p2_corr = flip(MTRasym_0p5uT_glu_pH7p2_corr);
fMTRasym_0p5uT_glu_pH7p0_corr = flip(MTRasym_0p5uT_glu_pH7p0_corr);
fMTRasym_0p5uT_glu_pH6p5_corr = flip(MTRasym_0p5uT_glu_pH6p5_corr);
fMTRasym_0p5uT_MnCl2_0p04mM_corr = flip(MTRasym_0p5uT_MnCl2_0p04mM_corr);
fMTRasym_0p5uT_MnCl2_0p08mM_corr = flip(MTRasym_0p5uT_MnCl2_0p08mM_corr);

%% Plot Z-spectra
figure(1)
hold on
plot(flip(d_Zspectra_0p5uT_eggwhite_corr))
plot(flip(d_Zspectra_0p5uT_glu_pH7p2_corr))
plot(flip(d_Zspectra_0p5uT_glu_pH7p0_corr))
plot(flip(d_Zspectra_0p5uT_glu_pH6p5_corr))
plot(flip(d_Zspectra_0p5uT_MnCl2_0p04mM_corr))
plot(flip(d_Zspectra_0p5uT_MnCl2_0p08mM_corr))

plot(flip(d_Zspectra_1p0uT_eggwhite_corr))
plot(flip(d_Zspectra_1p0uT_glu_pH7p2_corr))
plot(flip(d_Zspectra_1p0uT_glu_pH7p0_corr))
plot(flip(d_Zspectra_1p0uT_glu_pH6p5_corr))
plot(flip(d_Zspectra_1p0uT_MnCl2_0p04mM_corr))
plot(flip(d_Zspectra_1p0uT_MnCl2_0p08mM_corr))
xlabel('RF offset (ppm)');
ylabel('S/S_0 (%)')
ylim([0 1]);

%% Plot DSP Z-spectra
figure(2)
hold on
plot(flip(d_Zspectra_0p5uT_eggwhite_corr))
plot(flip(d_Zspectra_0p5uT_glu_pH7p2_corr))
plot(flip(d_Zspectra_0p5uT_glu_pH7p0_corr))
plot(flip(d_Zspectra_0p5uT_glu_pH6p5_corr))
plot(flip(d_Zspectra_0p5uT_MnCl2_0p04mM_corr))
plot(flip(d_Zspectra_0p5uT_MnCl2_0p08mM_corr))

plot(flip(Zspectra_DSP_eggwhite_0p5uT_corr))
plot(flip(Zspectra_DSP_glu_pH7p2_0p5uT_corr))
plot(flip(Zspectra_DSP_glu_pH7p0_0p5uT_corr))
plot(flip(Zspectra_DSP_glu_pH6p5_0p5uT_corr))
plot(flip(Zspectra_DSP_MnCl2_0p04mM_0p5uT_corr))
plot(flip(Zspectra_DSP_MnCl2_0p08mM_0p5uT_corr))
xlabel('RF offset (ppm)');
ylabel('S/S_0 (%)')
ylim([0.75 1]);

%% Plot MTR corrected Z-spectra
figure(3)
hold on
plot(flip(Zspectra_MTR_eggwhiteT_corr - Zspectra_MTR_eggwhiteT_corr(44)))
plot(flip(Zspectra_MTR_glu_pH7p2T_corr - Zspectra_MTR_glu_pH7p2T_corr(44)))
plot(flip(Zspectra_MTR_glu_pH7p0T_corr - Zspectra_MTR_glu_pH7p0T_corr(44)))
plot(flip(Zspectra_MTR_glu_pH6p5T_corr - Zspectra_MTR_glu_pH6p5T_corr(44)))
plot(flip(Zspectra_MTR_MnCl2_0p04mMT_corr - Zspectra_MTR_MnCl2_0p04mMT_corr(44)))
plot(flip(Zspectra_MTR_MnCl2_0p08mMT_corr - Zspectra_MTR_MnCl2_0p08mMT_corr(44)))

plot(fMTRasym_0p5uT_glu_pH7p2_corr)
plot(fMTRasym_0p5uT_glu_pH7p0_corr)
plot(fMTRasym_0p5uT_glu_pH6p5_corr)
plot(flip(1 - Zspectra_0p5uT_MnCl2_0p04mM_corr))
plot(flip(1 - Zspectra_0p5uT_MnCl2_0p08mM_corr))
xlabel('RF offset (ppm)');
ylabel('Corrected MTR_{DSP} (%)')
ylim([0 0.05]);

%% Plot AREX corrected Z-spectra
figure(4)
hold on
plot(flip(Zspectra_AREX_eggwhiteT_corr - Zspectra_AREX_eggwhiteT_corr(44)))
plot(flip(Zspectra_AREX_glu_pH7p2T_corr - Zspectra_AREX_glu_pH7p2T_corr(44)))
plot(flip(Zspectra_AREX_glu_pH7p0T_corr - Zspectra_AREX_glu_pH7p0T_corr(44)))
plot(flip(Zspectra_AREX_glu_pH6p5T_corr - Zspectra_AREX_glu_pH6p5T_corr(44)))
plot(flip(Zspectra_AREX_MnCl2_0p04mMT_corr - Zspectra_AREX_MnCl2_0p04mMT_corr(44)))
plot(flip(Zspectra_AREX_MnCl2_0p08mMT_corr - Zspectra_AREX_MnCl2_0p08mMT_corr(44)))
xlabel('RF offset (ppm)');
ylabel('Corrected AREX_{DSP} (%s^{-1})')
ylim([0 0.018]);
