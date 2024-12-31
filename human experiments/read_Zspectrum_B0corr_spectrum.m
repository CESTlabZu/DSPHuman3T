clear all; % Clear all variables from the workspace
close all; % Close all figure windows
clc;       % Clear the command window

% Matrix size
nx = 144; % Number of columns in the matrix
ny = 144; % Number of rows in the matrix

% Low saturation power and high saturation power
wl = 0.5; % Low saturation power
wh = 1;   % High saturation power

% Initialize matrices for storing results
Image_MTR_2pool(nx, ny, 47, 8) = 0;
Image_AREX_2pool(nx, ny, 47, 8) = 0;

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

% Loop through each subject
for sub_num = 1:6
    % Load the  brain data for the current subject
    load(sprintf('Sub%d/roi_brain.mat', sub_num))
    
    % generate WM/GM masks for the high resolution T2-weighted anatomy image
    roi_brain_anatomy_WM=double(imbinarize(imresize(roi_brain_WM, [512, 512])));
    roi_brain_anatomy_GM=double(imbinarize(imresize(roi_brain_GM, [512, 512])));

    % Calculate R1 map
    R1_map=1000./double(Imag1_T1);
    R1_map(R1_map==inf)=nan;
    R1_map(R1_map==-inf)=nan;
    Value_R1_WM(sub_num)=nansum(nansum(R1_map.*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Value_R1_GM(sub_num)=nansum(nansum(R1_map.*roi_brain_GM))./nansum(nansum(roi_brain_GM));


    % Read and normalize CEST data
    for ii = 1:51
        nImag1_0p5uT(:, :, ii) = double(Imag1_0p5uT(:, :, ii)) ./ mean(Imag1_0p5uT(:, :, [1, 2, 50, 51]), 3);
        nImag1_1uT(:, :, ii) = double(Imag1_1uT(:, :, ii)) ./ mean(Imag1_1uT(:, :, [1, 2, 50, 51]), 3);
    end

    % B0 correction for normalized images
    [nImag1_0p5uT_corr] = B0Correction(nImag1_0p5uT, roi_brain, FreqArray, range_cs, interp1NOE, deta_cs);
    [nImag1_1uT_corr] = B0Correction(nImag1_1uT, roi_brain, FreqArray, range_cs, interp1NOE, deta_cs);
    
    % Calculate DSP images
    nImag1_DSP_0p5uT_corr=1./(1+(1./nImag1_1uT_corr-1).*((wl)/(wh))^2);
    

    % LD mapping
    [Image_MTR_2pool, Image_AREX_2pool,  Image_reference_2pool] = ...
        LDmap(nImag1_0p5uT_corr, nx, ny, x_2pool, x, R1_map);

    nImag1_0p5uT_corr(nImag1_0p5uT_corr==0)=nan;
    nImag1_1uT_corr(nImag1_1uT_corr==0)=nan;
    nImag1_DSP_0p5uT_corr(nImag1_DSP_0p5uT_corr==0)=nan;
    Image_MTR_2pool(Image_MTR_2pool==0)=nan;
    Image_AREX_2pool(Image_AREX_2pool==0)=nan;    
    Image_reference_2pool(Image_reference_2pool==0)=nan;  
    T2_weighted_map(T2_weighted_map==0)=nan;  
    
    % Z-spectra with B0 correction
for ii=1:47
Zspectra_0p5uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_0p5uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
Zspectra_1uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_1uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

% DSP spectra with B0 correction
for ii=1:47
Zspectra_DSP_0p5uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_DSP_0p5uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_DSP_0p5uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_DSP_0p5uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end
for ii=1:47
Zspectra_MTR_WM_corr(ii,sub_num)=-nansum(nansum(squeeze(nImag1_0p5uT_corr(:,:,ii)-nImag1_DSP_0p5uT_corr(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_MTR_GM_corr(ii,sub_num)=-nansum(nansum(squeeze(nImag1_0p5uT_corr(:,:,ii)-nImag1_DSP_0p5uT_corr(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end
AREX_nImag1_0p5uT_corr_sub1=squeeze(1./nImag1_0p5uT_corr-1./nImag1_DSP_0p5uT_corr).*1000./double(Imag1_T1);
AREX_nImag1_0p5uT_corr_sub1(AREX_nImag1_0p5uT_corr_sub1==inf)=nan;
AREX_nImag1_0p5uT_corr_sub1(AREX_nImag1_0p5uT_corr_sub1==-inf)=nan;
for ii=1:47
Zspectra_AREX_WM_corr(ii,sub_num)=nansum(nansum(AREX_nImag1_0p5uT_corr_sub1(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_AREX_GM_corr(ii,sub_num)=nansum(nansum(AREX_nImag1_0p5uT_corr_sub1(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

% residual MT corrected DSP spectra with B0 correction
for ii=1:47
Zspectra_MTR_WM_corr_MT_corr(ii,sub_num)=-nansum(nansum(squeeze((nImag1_0p5uT_corr(:,:,ii)-nImag1_DSP_0p5uT_corr(:,:,ii))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_MTR_GM_corr_MT_corr(ii,sub_num)=-nansum(nansum(squeeze((nImag1_0p5uT_corr(:,:,ii)-nImag1_DSP_0p5uT_corr(:,:,ii))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end
AREX_nImag1_0p5uT_corr_MT_corr=squeeze((1./nImag1_0p5uT_corr(:,:,:)-1./nImag1_DSP_0p5uT_corr(:,:,:))-(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1);
AREX_nImag1_0p5uT_corr_MT_corr(AREX_nImag1_0p5uT_corr_MT_corr==inf)=nan;
AREX_nImag1_0p5uT_corr_MT_corr(AREX_nImag1_0p5uT_corr_MT_corr==-inf)=nan;
for ii=1:47
Zspectra_AREX_WM_corr_MT_corr(ii,sub_num)=nansum(nansum(AREX_nImag1_0p5uT_corr_MT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_AREX_GM_corr_MT_corr(ii,sub_num)=nansum(nansum(AREX_nImag1_0p5uT_corr_MT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

% asym spectra
for ii=1:23
MTRasym_0p5uT_WM_corr(ii,sub_num)=Zspectra_0p5uT_WM_corr(24-ii,sub_num)-Zspectra_0p5uT_WM_corr(24+ii,sub_num);
MTRasym_0p5uT_GM_corr(ii,sub_num)=Zspectra_0p5uT_GM_corr(24-ii,sub_num)-Zspectra_0p5uT_GM_corr(24+ii,sub_num);
end
fMTRasym_0p5uT_WM_corr=flip(MTRasym_0p5uT_WM_corr);
fMTRasym_0p5uT_GM_corr=flip(MTRasym_0p5uT_GM_corr);
fMTRasym_0p5uT_WM_corr(24,sub_num)=0;
fMTRasym_0p5uT_GM_corr(24,sub_num)=0;

for ii=1:23
AREXasym_0p5uT_WM_corr(ii,sub_num)=-(1./Zspectra_0p5uT_WM_corr(24-ii,sub_num)-1./Zspectra_0p5uT_WM_corr(24+ii,sub_num)).*Value_R1_WM(sub_num);
AREXasym_0p5uT_GM_corr(ii,sub_num)=-(1./Zspectra_0p5uT_GM_corr(24-ii,sub_num)-1./Zspectra_0p5uT_GM_corr(24+ii,sub_num)).*Value_R1_GM(sub_num);
end
fAREXasym_0p5uT_WM_corr=flip(AREXasym_0p5uT_WM_corr);
fAREXasym_0p5uT_GM_corr=flip(AREXasym_0p5uT_GM_corr);
fAREXasym_0p5uT_WM_corr(24,sub_num)=0;
fAREXasym_0p5uT_GM_corr(24,sub_num)=0;



for ii=1:47
Zspectra_MTR_2pool_WM(ii,sub_num)=nansum(nansum(Image_MTR_2pool(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_MTR_2pool_GM(ii,sub_num)=nansum(nansum(Image_MTR_2pool(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
Zspectra_AREX_2pool_WM(ii,sub_num)=nansum(nansum(Image_AREX_2pool(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_AREX_2pool_GM(ii,sub_num)=nansum(nansum(Image_AREX_2pool(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
Zspectra_reference_2pool_WM(ii,sub_num)=nansum(nansum(Image_reference_2pool(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_reference_2pool_GM(ii,sub_num)=nansum(nansum(Image_reference_2pool(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end



Value_WM_anatomy(sub_num)=nanmean(nanmean(nanmean(double(T2_weighted_map).*roi_brain_anatomy_WM)))./nanmean(nanmean(roi_brain_anatomy_WM));
Value_WM_MTR_corr_APT(sub_num)=squeeze(Zspectra_MTR_WM_corr(38, sub_num)-Zspectra_MTR_WM_corr(44, sub_num));
Value_WM_MTR_corr_NOE(sub_num)=squeeze(Zspectra_MTR_WM_corr(10, sub_num)-Zspectra_MTR_WM_corr(44, sub_num));
Value_WM_AREX_corr_APT(sub_num)=squeeze(Zspectra_AREX_WM_corr(38, sub_num)-Zspectra_AREX_WM_corr(44, sub_num));
Value_WM_AREX_corr_NOE(sub_num)=squeeze(Zspectra_AREX_WM_corr(10, sub_num)-Zspectra_AREX_WM_corr(44, sub_num));
Value_WM_MTR_2pool_APT(sub_num)=squeeze(Zspectra_MTR_2pool_WM(38, sub_num));
Value_WM_MTR_2pool_NOE(sub_num)=squeeze(Zspectra_MTR_2pool_WM(10, sub_num));
Value_WM_AREX_2pool_APT(sub_num)=squeeze(Zspectra_AREX_2pool_WM(38, sub_num));
Value_WM_AREX_2pool_NOE(sub_num)=squeeze(Zspectra_AREX_2pool_WM(10, sub_num));
Value_WM_MTRasym_corr(sub_num)=squeeze(fMTRasym_0p5uT_WM_corr(10, sub_num));
Value_WM_AREXasym_corr(sub_num)=squeeze(fAREXasym_0p5uT_WM_corr(10, sub_num));

Value_GM_anatomy(sub_num)=nanmean(nanmean(nanmean(double(T2_weighted_map).*roi_brain_anatomy_GM)))./nanmean(nanmean(roi_brain_anatomy_GM));
Value_GM_MTR_corr_APT(sub_num)=squeeze(Zspectra_MTR_GM_corr(38, sub_num)-Zspectra_MTR_GM_corr(44, sub_num));
Value_GM_MTR_corr_NOE(sub_num)=squeeze(Zspectra_MTR_GM_corr(10, sub_num)-Zspectra_MTR_GM_corr(44, sub_num));
Value_GM_AREX_corr_APT(sub_num)=squeeze(Zspectra_AREX_GM_corr(38, sub_num)-Zspectra_AREX_GM_corr(44, sub_num));
Value_GM_AREX_corr_NOE(sub_num)=squeeze(Zspectra_AREX_GM_corr(10, sub_num)-Zspectra_AREX_GM_corr(44, sub_num));
Value_GM_MTR_2pool_APT(sub_num)=squeeze(Zspectra_MTR_2pool_GM(38, sub_num));
Value_GM_MTR_2pool_NOE(sub_num)=squeeze(Zspectra_MTR_2pool_GM(10, sub_num));
Value_GM_AREX_2pool_APT(sub_num)=squeeze(Zspectra_AREX_2pool_GM(38, sub_num));
Value_GM_AREX_2pool_NOE(sub_num)=squeeze(Zspectra_AREX_2pool_GM(10, sub_num));
Value_GM_MTRasym_corr(sub_num)=squeeze(fMTRasym_0p5uT_GM_corr(10, sub_num));
Value_GM_AREXasym_corr(sub_num)=squeeze(fAREXasym_0p5uT_GM_corr(10, sub_num));



end





sub_range=[1:6]
% Fig. 6a
figure (1)
hold on
errorbar(flip(nanmean(Zspectra_0p5uT_WM_corr(:,sub_range)')), flip(nanstd(Zspectra_0p5uT_WM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_1uT_WM_corr(:,sub_range)')), flip(nanstd(Zspectra_1uT_WM_corr(:,sub_range)')));
ylim([0 1]);
xlabel('RF offset (ppm)');
ylabel('S/S_0 (%)');

% Fig. 6c
figure (2)
hold on
errorbar(flip(nanmean(Zspectra_0p5uT_WM_corr(:,sub_range)')), flip(nanstd(Zspectra_0p5uT_WM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_DSP_0p5uT_WM_corr(:,sub_range)')), flip(nanstd(Zspectra_DSP_0p5uT_WM_corr(:,sub_range)')));
ylim([0 1]);
xlabel('RF offset (ppm)');
ylabel('S/S_0 (%)');

% Fig. 6e
figure (3)
hold on
errorbar(flip(nanmean(Zspectra_MTR_WM_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_WM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_MTR_WM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_WM_corr_MT_corr(:,sub_range)')));
ylim([0 0.1]);
xlabel('RF offset (ppm)');
ylabel('MTR_{DSP} (%)');

% Fig. 6g
figure (4)
hold on
errorbar(flip(nanmean(Zspectra_AREX_WM_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_WM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_AREX_WM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_WM_corr_MT_corr(:,sub_range)')));
ylim([0 0.24]);
xlabel('RF offset (ppm)');
ylabel('AREX_{DSP} (%s^{-1})');

% Fig. 6b
figure (5)
hold on
errorbar(flip(nanmean(Zspectra_0p5uT_GM_corr(:,sub_range)')), flip(nanstd(Zspectra_0p5uT_GM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_1uT_GM_corr(:,sub_range)')), flip(nanstd(Zspectra_1uT_GM_corr(:,sub_range)')));
ylim([0 1]);
xlabel('RF offset (ppm)');
ylabel('S/S_0 (%)');

% Fig. 6d
figure (6)
hold on
errorbar(flip(nanmean(Zspectra_0p5uT_GM_corr(:,sub_range)')), flip(nanstd(Zspectra_0p5uT_GM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_DSP_0p5uT_GM_corr(:,sub_range)')), flip(nanstd(Zspectra_DSP_0p5uT_GM_corr(:,sub_range)')));
ylim([0 1]);
xlabel('RF offset (ppm)');
ylabel('S/S_0 (%)');

% Fig. 6f
figure (7)
hold on
errorbar(flip(nanmean(Zspectra_MTR_GM_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_GM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_MTR_GM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_GM_corr_MT_corr(:,sub_range)')));
ylim([0 0.1]);
xlabel('RF offset (ppm)');
ylabel('MTR_{DSP} (%)');

% Fig. 6h
figure (8)
hold on
errorbar(flip(nanmean(Zspectra_AREX_GM_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_GM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_AREX_GM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_GM_corr_MT_corr(:,sub_range)')));
ylim([0 0.24]);
xlabel('RF offset (ppm)');
ylabel('AREX_{DSP} (%s^{-1})');








% Fig. 7a
figure (9)
hold on
errorbar(flip(nanmean(Zspectra_0p5uT_WM_corr(:,sub_range)')), flip(nanstd(Zspectra_0p5uT_WM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_reference_2pool_WM(:,sub_range)')), flip(nanstd(Zspectra_reference_2pool_WM(:,sub_range)')));
ylim([0 1]);
xlabel('RF offset (ppm)');
ylabel('S/S_0 (%)');

% Fig. 7c
figure (10)
hold on
errorbar(flip(nanmean(Zspectra_MTR_2pool_WM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_WM(:,sub_range)')));
ylim([0 0.1]);
xlabel('RF offset (ppm)');
ylabel('MTR_{LD} (%)');

% Fig. 7e
figure (11)
hold on
errorbar(flip(nanmean(Zspectra_AREX_2pool_WM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_WM(:,sub_range)')));
ylim([0 0.16]);
xlabel('RF offset (ppm)');
ylabel('AREX_{LD} (%s^{-1})');

% Fig. 7g
figure (12)
hold on
errorbar((nanmean(fMTRasym_0p5uT_WM_corr(:,sub_range)')), (nanstd(fMTRasym_0p5uT_WM_corr(:,sub_range)')));
ylim([-0.08 0]);
xlabel('RF offset (ppm)');
ylabel('MTR_{LD} (%s)');


% Fig. 7i
figure (13)
hold on
errorbar((nanmean(fAREXasym_0p5uT_WM_corr(:,sub_range)')), (nanstd(fAREXasym_0p5uT_WM_corr(:,sub_range)')));
ylim([-0.15 0]);
xlabel('RF offset (ppm)');
ylabel('AREX_{LD} (%s^{-1})');



% Fig. 7b
figure (14)
hold on
errorbar(flip(nanmean(Zspectra_0p5uT_GM_corr(:,sub_range)')), flip(nanstd(Zspectra_0p5uT_GM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_reference_2pool_GM(:,sub_range)')), flip(nanstd(Zspectra_reference_2pool_GM(:,sub_range)')));
ylim([0 1]);
xlabel('RF offset (ppm)');
ylabel('S/S_0 (%)');

% Fig. 7d
figure (15)
hold on
errorbar(flip(nanmean(Zspectra_MTR_2pool_GM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_GM(:,sub_range)')));
ylim([0 0.1]);
xlabel('RF offset (ppm)');
ylabel('MTR_{LD} (%)');

% Fig. 7f
figure (16)
hold on
errorbar(flip(nanmean(Zspectra_AREX_2pool_GM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_GM(:,sub_range)')));
ylim([0 0.16]);
xlabel('RF offset (ppm)');
ylabel('AREX_{LD} (%s^{-1})');


% Fig. 7h
figure (17)
hold on
errorbar((nanmean(fMTRasym_0p5uT_GM_corr(:,sub_range)')), (nanstd(fMTRasym_0p5uT_GM_corr(:,sub_range)')));
ylim([-0.08 0]);
xlabel('RF offset (ppm)');
ylabel('MTR_{LD} (%s)');

% Fig. 7j
figure (18)
hold on
errorbar((nanmean(fAREXasym_0p5uT_GM_corr(:,sub_range)')), (nanstd(fAREXasym_0p5uT_GM_corr(:,sub_range)')));
ylim([-0.15 0]);
xlabel('RF offset (ppm)');
ylabel('AREX_{LD} (%s^{-1})');




%Fig.8a
figure (19)
hold on
errorbar(flip(nanmean(Zspectra_MTR_WM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_WM_corr_MT_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_MTR_2pool_WM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_WM(:,sub_range)')));
errorbar((nanmean(fMTRasym_0p5uT_WM_corr(:,sub_range)')), (nanstd(fMTRasym_0p5uT_WM_corr(:,sub_range)')));
ylim([-0.06 0.07]);
xlabel('RF offset (ppm)');
ylabel('MTR (%)');

%Fig.8b
figure (20)
hold on
errorbar(flip(nanmean(Zspectra_MTR_GM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_GM_corr_MT_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_MTR_2pool_GM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_GM(:,sub_range)')));
errorbar((nanmean(fMTRasym_0p5uT_GM_corr(:,sub_range)')), (nanstd(fMTRasym_0p5uT_GM_corr(:,sub_range)')));
ylim([-0.06 0.07]);
xlabel('RF offset (ppm)');
ylabel('MTR (%)');

%Fig.8c
figure (21)
hold on
errorbar(flip(nanmean(Zspectra_AREX_WM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_WM_corr_MT_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_AREX_2pool_WM(:,sub_range)')), flip(nanstd(Zspectra_AREX_2pool_WM(:,sub_range)')));
errorbar((nanmean(fAREXasym_0p5uT_WM_corr(:,sub_range)')), (nanstd(fAREXasym_0p5uT_WM_corr(:,sub_range)')));
ylim([-0.15 0.15]);
xlabel('RF offset (ppm)');
ylabel('AREX (%s^{-1})');

%Fig.8d
figure (22)
hold on
errorbar(flip(nanmean(Zspectra_AREX_GM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_GM_corr_MT_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_AREX_2pool_GM(:,sub_range)')), flip(nanstd(Zspectra_AREX_2pool_GM(:,sub_range)')));
errorbar((nanmean(fAREXasym_0p5uT_GM_corr(:,sub_range)')), (nanstd(fAREXasym_0p5uT_GM_corr(:,sub_range)')));
ylim([-0.15 0.15]);
xlabel('RF offset (ppm)');
ylabel('AREX (%s^{-1})');






% Fig.10a
% statistical results   anatomy
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_anatomy(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_anatomy(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(23);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_anatomy(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_anatomy(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(23);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_GM_anatomy(sub_range), Value_WM_anatomy(sub_range))
ylabel('T_2 weighted anatomy (a.u.)')
ylim([0 800]);
path = sprintf('%s%d.fig','FigS7\',1);
set(gca,'xtick',[]);set(gca,'ytick',[]);





% Fig.10b
% statistical results   R1
x_zhu=[1,2]
y_zhu(1)=mean(Value_R1_WM(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_R1_WM(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(24);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_R1_GM(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_R1_GM(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(24);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_R1_WM(sub_range), Value_R1_GM(sub_range))
ylabel('R_{1obs} (s^{-1})')
ylim([0 2]);
path = sprintf('%s%d.fig','FigS7\',2);
set(gca,'xtick',[]);set(gca,'ytick',[]);






% Fig.10c
% statistical results   MTRdsp APT
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_MTR_corr_APT(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_MTR_corr_APT(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(25);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_MTR_corr_APT(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_MTR_corr_APT(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(25);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_GM_MTR_corr_APT(sub_range), Value_WM_MTR_corr_APT(sub_range))
ylabel('MTR_{DSP} at 3.5ppm (%)')
ylim([0 0.04]);
path = sprintf('%s%d.fig','FigS7\',3);
set(gca,'xtick',[]);set(gca,'ytick',[]);






% Fig.10d
% statistical results   AREXdsp APT
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_AREX_corr_APT(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_AREX_corr_APT(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(26);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_AREX_corr_APT(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_AREX_corr_APT(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(26);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_AREX_corr_APT(sub_range), Value_GM_AREX_corr_APT(sub_range))
ylabel('AREX_{DSP} at 3.5ppm (%s^{-1})')
ylim([0 0.06]);
path = sprintf('%s%d.fig','FigS7\',4);
set(gca,'xtick',[]);set(gca,'ytick',[]);










% Fig.10e
% statistical results   MTRld APT
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_MTR_2pool_APT(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_MTR_2pool_APT(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(27);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_MTR_2pool_APT(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_MTR_2pool_APT(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(27);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_MTR_2pool_APT(sub_range), Value_GM_MTR_2pool_APT(sub_range))
ylabel('MTR_{LD} at 3.5ppm (%)')
ylim([0 0.08]);
path = sprintf('%s%d.fig','FigS7\',5);
set(gca,'xtick',[]);set(gca,'ytick',[]);








% Fig.10f
% statistical results   AREXld APT
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_AREX_2pool_APT(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_AREX_2pool_APT(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(28);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_AREX_2pool_APT(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_AREX_2pool_APT(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(28);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_AREX_2pool_APT(sub_range), Value_GM_AREX_2pool_APT(sub_range))
ylabel('AREX_{LD} at 3.5ppm (%s^{-1})')
ylim([0 0.12]);
path = sprintf('%s%d.fig','FigS7\',6);
set(gca,'xtick',[]);set(gca,'ytick',[]);






% Fig.10g
% statistical results   MTRasym
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_MTRasym_corr(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_MTRasym_corr(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(29);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_MTRasym_corr(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_MTRasym_corr(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(29);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_MTRasym_corr(sub_range), Value_GM_MTRasym_corr(sub_range))
ylabel('MTR_{asym} at 3.5ppm (%)')
ylim([-0.1 0]);
path = sprintf('%s%d.fig','FigS7\',7);
set(gca,'xtick',[]);set(gca,'ytick',[]);






% Fig.10h
% statistical results   AREXasym
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_AREXasym_corr(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_AREXasym_corr(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(30);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_AREXasym_corr(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_AREXasym_corr(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(30);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_AREXasym_corr(sub_range), Value_GM_AREXasym_corr(sub_range))
ylabel('AREX_{asym} at 3.5ppm (%s^{-1})')
ylim([-0.2 0]);
path = sprintf('%s%d.fig','FigS7\',8);
set(gca,'xtick',[]);set(gca,'ytick',[]);






% Fig.10i
% statistical results   MTRdsp NOE
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_MTR_corr_NOE(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_MTR_corr_NOE(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(31);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_MTR_corr_NOE(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_MTR_corr_NOE(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(31);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_MTR_corr_NOE(sub_range), Value_GM_MTR_corr_NOE(sub_range))
ylabel('MTR_{DSP} at -3.5ppm (%)')
ylim([0 0.1]);
path = sprintf('%s%d.fig','FigS7\',9);
set(gca,'xtick',[]);set(gca,'ytick',[]);







% Fig.10j
% statistical results   AREXdsp NOE
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_AREX_corr_NOE(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_AREX_corr_NOE(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(32);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_AREX_corr_NOE(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_AREX_corr_NOE(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(32);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_AREX_corr_NOE(sub_range), Value_GM_AREX_corr_NOE(sub_range))
ylabel('AREX_{DSP} at -3.5ppm (%s^{-1})')
ylim([0 0.15]);
path = sprintf('%s%d.fig','FigS7\',10);
set(gca,'xtick',[]);set(gca,'ytick',[]);







% Fig.10k
% statistical results   MTRld NOE
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_MTR_2pool_NOE(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_MTR_2pool_NOE(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(33);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_MTR_2pool_NOE(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_MTR_2pool_NOE(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(33);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_MTR_2pool_NOE(sub_range), Value_GM_MTR_2pool_NOE(sub_range))
ylabel('MTR_{LD} at -3.5ppm (%)')
ylim([0 0.1]);
path = sprintf('%s%d.fig','FigS7\',11);
set(gca,'xtick',[]);set(gca,'ytick',[]);








% Fig.10l
% statistical results   AREXld NOE
x_zhu=[1,2]
y_zhu(1)=mean(Value_WM_AREX_2pool_NOE(sub_range))
y_zhu(2)=0;
bb_zhu(1)=std(Value_WM_AREX_2pool_NOE(sub_range))
bb_zhu(2)=0;

lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;
L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(34);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

x_zhu=[1,2]
y_zhu(1)=0;
y_zhu(2)=mean(Value_GM_AREX_2pool_NOE(sub_range))
bb_zhu(1)=0;
bb_zhu(2)=std(Value_GM_AREX_2pool_NOE(sub_range))
lower_zhu = y_zhu - bb_zhu;
upper_zhu = y_zhu + bb_zhu;

L_zhu = y_zhu - lower_zhu;
U_zhu = upper_zhu -y_zhu;

figure(34);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_AREX_2pool_APT(sub_range), Value_GM_AREX_2pool_APT(sub_range))
ylabel('AREX_{LD} at -3.5ppm (%s^{-1})')
ylim([0 0.15]);
path = sprintf('%s%d.fig','FigS7\',12);
set(gca,'xtick',[]);set(gca,'ytick',[]);







% Fig. S27
Value_GMWM_AREX_corr_APT=[Value_GM_AREX_corr_APT Value_WM_AREX_corr_APT];
Value_GMWM_AREX_corr_NOE=[Value_GM_AREX_corr_NOE Value_WM_AREX_corr_NOE];
Value_GMWM_MTR_corr_APT=[Value_GM_MTR_corr_APT Value_WM_MTR_corr_APT];
Value_GMWM_MTR_corr_NOE=[Value_GM_MTR_corr_NOE Value_WM_MTR_corr_NOE];

Value_GMWM_AREX_2pool_APT=[Value_GM_AREX_2pool_APT Value_WM_AREX_2pool_APT];
Value_GMWM_AREX_2pool_NOE=[Value_GM_AREX_2pool_NOE Value_WM_AREX_2pool_NOE];
Value_GMWM_MTR_2pool_APT=[Value_GM_MTR_2pool_APT Value_WM_MTR_2pool_APT];
Value_GMWM_MTR_2pool_NOE=[Value_GM_MTR_2pool_NOE Value_WM_MTR_2pool_NOE];

Value_GMWM_AREXasym_corr=[Value_GM_AREXasym_corr Value_WM_AREXasym_corr];
Value_GMWM_MTRasym_corr=[Value_GM_MTRasym_corr Value_WM_MTRasym_corr];


figure (35) 
x=Value_GMWM_MTR_corr_APT*100;
y=Value_GMWM_MTR_2pool_APT*100;
x_GM=Value_GM_MTR_corr_APT*100;
y_GM=Value_GM_MTR_2pool_APT*100;
x_WM=Value_WM_MTR_corr_APT*100;
y_WM=Value_WM_MTR_2pool_APT*100;

scatter(x, y)
xlabel('MTR_{DSP} at 3.5ppm (%)');
ylabel('MTR_{LD} at 3.5ppm (%)');
xlim([0.5 3]);
ylim([2 5]);

p = polyfit(x, y, 1);
y_fit = polyval(p, x);
hold on;
plot(x, y_fit, 'k--', 'LineWidth', 2); % Plot the linear fit

scatter(x_GM, y_GM, 'b')
scatter(x_WM, y_WM, 'r')

[R, p_value] = corr(x', y', 'Type', 'Spearman');
correlation_coefficient = R(1, 1); % Extract the correlation coefficient
% Get the current axes limits
ax = gca;
x_limits = ax.XLim;
y_limits = ax.YLim;

% Set text position at the top-left corner
text_position_x = x_limits(1) + 0.05 * (x_limits(2) - x_limits(1)); % x position
text_position_y = y_limits(2) - 0.05 * (y_limits(2) - y_limits(1)); % y position
text_str = ['r = ', num2str(correlation_coefficient, '%.2f'), ', p = ', num2str(p_value, '%.2f')];

% Display the text
text(text_position_x, text_position_y, text_str, 'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w', 'VerticalAlignment', 'top');





% Figure supporting
figure (36) 
x=Value_GMWM_AREX_corr_APT*100;
y=Value_GMWM_AREX_2pool_APT*100;
x_GM=Value_GM_AREX_corr_APT*100;
y_GM=Value_GM_AREX_2pool_APT*100;
x_WM=Value_WM_AREX_corr_APT*100;
y_WM=Value_WM_AREX_2pool_APT*100;

scatter(x, y)
xlabel('AREX_{DSP} at 3.5ppm (%)');
ylabel('AREX_{LD} at 3.5ppm (%)');
xlim([1 5]);
ylim([2 8]);

p = polyfit(x, y, 1);
y_fit = polyval(p, x);
hold on;
plot(x, y_fit, 'k--', 'LineWidth', 2); % Plot the linear fit

scatter(x_GM, y_GM, 'b')
scatter(x_WM, y_WM, 'r')

[R, p_value] = corr(x', y', 'Type', 'Spearman');
correlation_coefficient = R(1, 1); % Extract the correlation coefficient
% Get the current axes limits
ax = gca;
x_limits = ax.XLim;
y_limits = ax.YLim;

% Set text position at the top-left corner
text_position_x = x_limits(1) + 0.05 * (x_limits(2) - x_limits(1)); % x position
text_position_y = y_limits(2) - 0.05 * (y_limits(2) - y_limits(1)); % y position
text_str = ['r = ', num2str(correlation_coefficient, '%.2f'), ', p = ', num2str(p_value, '%.2f')];

% Display the text
text(text_position_x, text_position_y, text_str, 'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w', 'VerticalAlignment', 'top');







% Figure supporting
figure (37) 
x=Value_GMWM_MTR_corr_NOE*100;
y=Value_GMWM_MTR_2pool_NOE*100;
x_GM=Value_GM_MTR_corr_NOE*100;
y_GM=Value_GM_MTR_2pool_NOE*100;
x_WM=Value_WM_MTR_corr_NOE*100;
y_WM=Value_WM_MTR_2pool_NOE*100;

scatter(x, y)
xlabel('MTR_{DSP} at -3.5ppm (%)');
ylabel('MTR_{LD} at -3.5ppm (%)');
xlim([3 6]);
ylim([4 8]);

p = polyfit(x, y, 1);
y_fit = polyval(p, x);
hold on;
plot(x, y_fit, 'k--', 'LineWidth', 2); % Plot the linear fit

scatter(x_GM, y_GM, 'b')
scatter(x_WM, y_WM, 'r')

[R, p_value] = corr(x', y', 'Type', 'Spearman');
correlation_coefficient = R(1, 1); % Extract the correlation coefficient
% Get the current axes limits
ax = gca;
x_limits = ax.XLim;
y_limits = ax.YLim;

% Set text position at the top-left corner
text_position_x = x_limits(1) + 0.05 * (x_limits(2) - x_limits(1)); % x position
text_position_y = y_limits(2) - 0.05 * (y_limits(2) - y_limits(1)); % y position
text_str = ['r = ', num2str(correlation_coefficient, '%.2f'), ', p = ', num2str(p_value, '%.2f')];

% Display the text
text(text_position_x, text_position_y, text_str, 'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w', 'VerticalAlignment', 'top');






% Figure supporting
figure (38) 
x=Value_GMWM_AREX_corr_NOE*100;
y=Value_GMWM_AREX_2pool_NOE*100;
x_GM=Value_GM_AREX_corr_NOE*100;
y_GM=Value_GM_AREX_2pool_NOE*100;
x_WM=Value_WM_AREX_corr_NOE*100;
y_WM=Value_WM_AREX_2pool_NOE*100;

scatter(x, y)
xlabel('AREX_{DSP} at -3.5ppm (%)');
ylabel('AREX_{LD} at -3.5ppm (%)');
xlim([2 14]);
ylim([4 16]);

p = polyfit(x, y, 1);
y_fit = polyval(p, x);
hold on;
plot(x, y_fit, 'k--', 'LineWidth', 2); % Plot the linear fit

scatter(x_GM, y_GM, 'b')
scatter(x_WM, y_WM, 'r')

[R, p_value] = corr(x', y', 'Type', 'Spearman');
correlation_coefficient = R(1, 1); % Extract the correlation coefficient
% Get the current axes limits
ax = gca;
x_limits = ax.XLim;
y_limits = ax.YLim;

% Set text position at the top-left corner
text_position_x = x_limits(1) + 0.05 * (x_limits(2) - x_limits(1)); % x position
text_position_y = y_limits(2) - 0.05 * (y_limits(2) - y_limits(1)); % y position
text_str = ['r = ', num2str(correlation_coefficient, '%.2f'), ', p = ', num2str(p_value, '%.2f')];

% Display the text
text(text_position_x, text_position_y, text_str, 'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w', 'VerticalAlignment', 'top');






