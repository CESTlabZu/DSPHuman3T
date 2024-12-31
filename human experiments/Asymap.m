function [Image_MTRasym_corr, Image_AREXasym_corr] = Asymap(nImag_corr, R1_map)
    % Function to perform asymmetry calculation

    % Calculate MTR asymmetry
    Image_MTRasym_corr = nImag_corr(:, :, 10) - nImag_corr(:, :, 38);
    
    % Calculate AREX asymmetry
    Image_AREXasym_corr = -(1 ./ nImag_corr(:, :, 10) - 1 ./ nImag_corr(:, :, 38)) .* R1_map;
end