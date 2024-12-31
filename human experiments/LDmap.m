function [Image_MTR_2pool, Image_AREX_2pool, Image_reference_2pool] = LDmap(nImag_corr, nx, ny, x_2pool, x, R1_map)
    % Function to perform Lorentzian difference calculation 

Image_MTR_2pool=zeros(nx,ny,length(x));
Image_AREX_2pool=zeros(nx,ny,length(x));
Image_reference_2pool=zeros(nx,ny,length(x));

    for xx = 1:nx
        for yy = 1:ny
            Zspectra_0p5uT_corr_single = squeeze(nImag_corr(xx, yy, :));

            if Zspectra_0p5uT_corr_single(1) == 0
                continue;
            end

            sig_2pool = zeros(1, 9);
            for i = 1:2
                sig_2pool(i) = 1 - Zspectra_0p5uT_corr_single(i);
            end
            for i = 3:7
                sig_2pool(i) = 1 - Zspectra_0p5uT_corr_single(i + 19);
            end
            for i = 8:9
                sig_2pool(i) = 1 - Zspectra_0p5uT_corr_single(i + 38);
            end


            beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
            lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
            ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

            Delta = [1]; % constants

            options = optimset('lsqcurvefit');
            options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8, 'MaxFunEvals', 5e4 * length(x), 'MaxIter', 2e5);

            [beta_2pool, ~, ~, ~, ~, ~, ~] = ...
                lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta);

            sig_simu_2pool = matsolv_2pool(beta_2pool, x, Delta);

            Image_MTR_2pool(xx, yy, :) = 1 - sig_simu_2pool - squeeze(Zspectra_0p5uT_corr_single);
            
               Image_AREX_2pool(xx, yy, :) = ...
                (1 ./ squeeze(Zspectra_0p5uT_corr_single) - 1 ./ (1 - sig_simu_2pool)) .* R1_map(xx, yy);
            
            Image_reference_2pool(xx, yy, :) = 1-sig_simu_2pool;
        end
    end
end