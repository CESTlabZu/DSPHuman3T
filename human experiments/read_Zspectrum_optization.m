clear all; % Clear all variables from the workspace
close all; % Close all figure windows
clc;       % Clear the command window

% Matrix size
nx = 144; ny = 144;

% Powers combinations for each subject
powers = {{'0p3uT', '0p6uT'}, {'0p4uT', '0p8uT'}, {'0p5uT', '1uT'}, {'0p6uT', '1p2uT'}, {'0p7uT', '1p4uT'}};
powers_value = {{0.3, 0.6}, {0.4, 0.8}, {0.5, 1.0}, {0.6, 1.2}, {0.7, 1.4}};

% ROI types
roi_types = {'WM', 'GM'};

% Normalize images
normalize_image = @(field) arrayfun(@(ii) double(field(:,:,ii)) ./ mean(field(:,:, [1,2,50,51]), 3), 1:51, 'UniformOutput', false);

% Calculate Z-spectra
calc_zspectra = @(fields, roi) arrayfun(@(ii) nansum(nansum(double(fields{ii}) .* roi)) / nansum(nansum(roi)), 1:51);

% Calculate DSP spectra
calc_zspectra_DSP = @(zspectra_low, zspectra_high, wl, wh) 1 ./ (1 + (1 ./ zspectra_high - 1) .* (wl / wh)^2);

% Initialize MTRDSP matrix
MTRDSP_matrix = zeros(6, 2, length(powers), 51);

% Loop through each subject
for sub_num = 1:6
    % Load the brain data for the current subject
    load(sprintf('Sub%d/roi_brain.mat', sub_num))
    
    for roi_type = 1:2 
        roi = eval(['roi_brain_' roi_types{roi_type}]);
        for power = 1:length(powers)
            wl = powers_value{power}{1};
            wh = powers_value{power}{2};   
            
            imag_low = eval(sprintf('Imag1_%s', powers{power}{1}));
            imag_high = eval(sprintf('Imag1_%s', powers{power}{2}));
    
            norm_low = normalize_image(imag_low);
            norm_high = normalize_image(imag_high);
         
            zspectra_low = calc_zspectra(norm_low, roi);
            zspectra_high = calc_zspectra(norm_high, roi);
            
            % Calculate DSP-corrected Z-spectra
            Zspectra_DSP = calc_zspectra_DSP(zspectra_low, zspectra_high, wl, wh);
            
            % Calculate MTRDSP
            MTRDSP = Zspectra_DSP - zspectra_low;
            
            % Calculate residual MT corrected MTRDSP
            MTRDSP_matrix(sub_num, roi_type, power, :) = MTRDSP - MTRDSP(46);
        end
    end
end

% Draw figure Fig. S19a
figure(1);
plot(flip(squeeze(nanmean(MTRDSP_matrix(:, 1, :, :)))'), 'LineWidth', 2);
xlim([6 46]);
ylim([0 0.05]);
xlabel('RF offset (ppm)');
ylabel('MTR_{DSP} (%)');
title('Mean MTR_{DSP} for WM');


% Draw figure Fig. S19b
figure(2);
plot(flip(squeeze(nanmean(MTRDSP_matrix(:, 2, :, :)))'), 'LineWidth', 2);
xlim([6 46]);
ylim([0 0.05]);
xlabel('RF offset (ppm)');
ylabel('MTR_{DSP} (%)');
title('Mean MTR_{DSP} for GM');

