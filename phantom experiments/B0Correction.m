function [nImag_corr] = B0Correction(nImag, roi_brain, FreqArray, range_cs, interp1NOE, deta_cs)
    % Function to perform B0 correction

    [Length, Width, Freq] = size(nImag);
    nImag_corr = zeros(Length, Width, 47);
    InterplCESTData = zeros(Length, Width, interp1NOE);

    % Initial Mask and Interpolation
    Mask1 = zeros(Length, Width, Freq);
    for xx = 1:Length
        for yy = 1:Width
            Raw_Z = squeeze(nImag(xx, yy, 3:49) .* roi_brain(xx, yy))';
            if Raw_Z(1) == 0 || isnan(Raw_Z(1))
                Mask1(xx, yy, :) = 0;
            else
                Mask1(xx, yy, :) = 1;
                AA = interp1(FreqArray, Raw_Z, range_cs, 'spline');
                InterplCESTData(xx, yy, :) = AA;
            end
        end
    end

    % Initialize B0shift map and Corrected Z
    B0shiftmap = zeros(Length, Width);
    Corrected_Z = zeros(Length, Width, interp1NOE);
    B0shiftAcquisitionNumber = interp1NOE;
    detacenter = 0;

    % Perform B0 Correction
    for xx = 1:Length
        for yy = 1:Width
            B0shiftvalue = squeeze(InterplCESTData(xx, yy, :));
            [minvalue, position] = min(B0shiftvalue);
            if position > 1 && position < B0shiftAcquisitionNumber
                deta = B0shiftvalue(position-1) - B0shiftvalue(position+1);
                if deta > 0
                    detacenter = 0.5 * deta_cs - (0.5 * deta_cs) / ((B0shiftvalue(position-1) - minvalue) / (B0shiftvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5 * deta_cs + (0.5 * deta_cs) / ((B0shiftvalue(position+1) - minvalue) / (B0shiftvalue(position-1) - minvalue));
                else
                    detacenter = 0;
                end
            end
            center = ((B0shiftAcquisitionNumber + 1) / 2 - position) * deta_cs - detacenter;
            B0shiftmap(xx, yy) = center;
            detacenter = 0;

            corrected_cs = range_cs + double(center);
            Uncorrected_Z = squeeze(InterplCESTData(xx, yy, :));
            if Uncorrected_Z(1) > 0.05 && abs(center) < 2
                Corrected_Z(xx, yy, :) = interp1(corrected_cs, Uncorrected_Z, range_cs, 'spline');
                TarValue_2 = interp1(corrected_cs, Uncorrected_Z, FreqArray, 'spline');
                nImag_corr(xx, yy, :) = TarValue_2;
            end

   end
    end
end
