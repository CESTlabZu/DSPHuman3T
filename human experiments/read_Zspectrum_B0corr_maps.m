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
    
    % Calculate R1 map
    R1_map = 1000 ./ double(Imag1_T1);
    R1_map(R1_map == inf) = nan;
    R1_map(R1_map == -inf) = nan;

    % Initialize corrected MTR and AREX images
    MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47) = 0;
    AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47) = 0;

    % Read and normalize CEST data
    for ii = 1:51
        nImag1_0p5uT(:, :, ii) = double(Imag1_0p5uT(:, :, ii)) ./ mean(Imag1_0p5uT(:, :, [1, 2, 50, 51]), 3);
        nImag1_1uT(:, :, ii) = double(Imag1_1uT(:, :, ii)) ./ mean(Imag1_1uT(:, :, [1, 2, 50, 51]), 3);
    end

    % B0 correction for normalized images
    [nImag1_0p5uT_corr] = B0Correction(nImag1_0p5uT, roi_brain, FreqArray, range_cs, interp1NOE, deta_cs);
    [nImag1_1uT_corr] = B0Correction(nImag1_1uT, roi_brain, FreqArray, range_cs, interp1NOE, deta_cs);

    % Residual MT corrected DSP maps
    [Image_MTR_APT_0p5uT_corr(:, :, sub_num), Image_MTR_NOE_0p5uT_corr(:, :, sub_num), ...
        Image_AREX_APT_0p5uT_corr(:, :, sub_num), Image_AREX_NOE_0p5uT_corr(:, :, sub_num)] = ...
        DSPMap(nImag1_0p5uT_corr, nImag1_1uT_corr, wl, wh, roi_brain, Imag1_T1);

    % LD mapping
    [Image_MTR_2pool(:, :, :, sub_num), Image_AREX_2pool(:, :, :, sub_num)] = ...
        LDmap(nImag1_0p5uT_corr, nx, ny, x_2pool, x, R1_map);

    % Asymmetry mapping
    [Image_MTRasym_0p5uT_corr(:, :, sub_num), Image_AREXasym_0p5uT_corr(:, :, sub_num)] = ...
        Asymap(nImag1_0p5uT_corr, R1_map);
end





% Draw images, change sub_num to select different subject
% #1-#6 in this code is #6-#1 in the paper 
load('Sub6\roi_brain')
sub_num=6
%Fig. 9a T2-weighted anatomy
figure (1)
imagesc(T2_weighted_map, [0 1000])
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'a', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('T_2-weighted anatomy (a.u.)', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',1);
savefig(path) 

%Fig. 9b T1
figure (2)
imagesc(1000./double(Imag1_T1).*roi_brain,[0 2])
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'b', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('R_{1obs} (s^{-1})', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',2);
savefig(path) 

% Fig. 9c MTRDSP at 3.5ppm
figure (3)
imagesc(Image_MTR_APT_0p5uT_corr(:,:,sub_num)*100,[0 0.04]*100)
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'c', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('MTR_{DSP} at 3.5ppm (%)', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',3);
savefig(path) 

% Fig. 9d AREXDSP at 3.5ppm
figure (4)
imagesc(Image_AREX_APT_0p5uT_corr(:,:,sub_num)*100,[0 0.06]*100)
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'd', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('AREX_{DSP} at 3.5ppm (%s^{-1})', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',4);
savefig(path) 

% Fig. 9e MTRLD at 3.5ppm
figure (5)
imagesc(Image_MTR_2pool(:,:,38,sub_num).*roi_brain*100, [0 0.08]*100)
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'e', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('MTR_{LD} at 3.5ppm (%)', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',5);
savefig(path) 

% Fig. 9f AREXLD at 3.5ppm
figure (6)
imagesc(Image_AREX_2pool(:,:,38,sub_num).*roi_brain*100, [0 0.12]*100)
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'f', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('AREX_{LD} at 3.5ppm (%s^{-1})', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',6);
savefig(path) 

% Fig. 9g MTRasym at 3.5ppm 
figure (7)
Image_MTRasym_0p5uT_corr(Image_MTRasym_0p5uT_corr==0)=nan;
imagesc(Image_MTRasym_0p5uT_corr(:,:,sub_num).*roi_brain*100,[-0.1 0]*100)
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'g', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('MTR_{asym} at 3.5ppm (%)', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',7);
savefig(path) 

% Fig. 9h AREXasym at 3.5ppm
figure (8)
Image_AREXasym_0p5uT_corr(Image_MTRasym_0p5uT_corr==0)=nan;
imagesc(Image_AREXasym_0p5uT_corr(:,:,sub_num).*roi_brain*100,[-0.2 0]*100)
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'h', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('AREX_{asym} at 3.5ppm (%s^{-1})', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',8);
savefig(path) 




% NOE
% Fig. 9i MTRDSP at -3.5ppm
figure (9)
imagesc(Image_MTR_NOE_0p5uT_corr(:,:,sub_num)*100,[0 0.1]*100)
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'i', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('MTR_{DSP} at -3.5ppm (%)', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',9);
savefig(path) 

% Fig. 9j MTRDSP at -3.5ppm
figure (10)
imagesc(Image_AREX_NOE_0p5uT_corr(:,:,sub_num)*100,[0 0.15]*100)
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'j', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('AREX_{DSP} at -3.5ppm (%s^{-1})', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',10);
savefig(path) 

% Fig. 9k MTRLD at -3.5ppm
figure (11)
imagesc(Image_MTR_2pool(:,:,10,sub_num).*roi_brain*100, [0 0.1]*100)
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'k', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('MTR_{LD} at -3.5ppm (%)', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',11);
savefig(path) 

% Fig. 9l AREXLD at -3.5ppm
figure (12)
imagesc(Image_AREX_2pool(:,:,10,sub_num).*roi_brain*100, [0 0.15]*100)
pos = get(gca, 'Position');
annotation('textbox', [pos(1) - 0.08, pos(2) - 0.0, 0.03, 0.03], 'String', 'l', ...
           'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none');
xlabel('AREX_{LD} at -3.5ppm (%s^{-1})', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
path = sprintf('%s%d.fig','figures\',12);
savefig(path) 


%Fig. S2
% mask WM
figure (13)
imagesc(roi_brain_WM)
colormap(gray);
xlabel('WM mask', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
% mask GM
figure (14)
imagesc(roi_brain_GM)
colormap(gray);
xlabel('GM mask', 'Interpreter', 'tex');
set(gca,'xtick',[]);set(gca,'ytick',[]);
