function [Image_MTR_APT_corr, Image_MTR_NOE_corr, Image_AREX_APT_corr, Image_AREX_NOE_corr] = DSPMap(nImag1_low_corr, nImag1_high_corr, wl, wh, roi_brain, Imag1_T1)
    % Function to calculate residual MT corrected DSP maps 

    % Calculate DSP corrected Z-spectra
    nImag1_DSP_low_corr = 1 ./ (1 + (1 ./ nImag1_high_corr - 1) .* ((wl / wh) ^ 2));

    % Calculate MTR APT and NOE corrected images
    Image_MTR_APT_corr = -(squeeze(nImag1_low_corr(:, :, 38) - nImag1_DSP_low_corr(:, :, 38)) - ...
                           squeeze(nImag1_low_corr(:, :, 44) - nImag1_DSP_low_corr(:, :, 44))) .* roi_brain;
                       
    Image_MTR_NOE_corr = -(squeeze(nImag1_low_corr(:, :, 10) - nImag1_DSP_low_corr(:, :, 10)) - ...
                           squeeze(nImag1_low_corr(:, :, 44) - nImag1_DSP_low_corr(:, :, 44))) .* roi_brain;

    % Calculate AREX APT and NOE corrected images
    Image_AREX_APT_corr = (squeeze(1 ./ nImag1_low_corr(:, :, 38) - 1 ./ nImag1_DSP_low_corr(:, :, 38)) - ...
                           squeeze(1 ./ nImag1_low_corr(:, :, 44) - 1 ./ nImag1_DSP_low_corr(:, :, 44))) .* 1000 ./ double(Imag1_T1) .* roi_brain;
                       
    Image_AREX_NOE_corr = (squeeze(1 ./ nImag1_low_corr(:, :, 10) - 1 ./ nImag1_DSP_low_corr(:, :, 10)) - ...
                           squeeze(1 ./ nImag1_low_corr(:, :, 44) - 1 ./ nImag1_DSP_low_corr(:, :, 44))) .* 1000 ./ double(Imag1_T1) .* roi_brain;            
end
