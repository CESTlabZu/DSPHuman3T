clear all; close all; clc;
% matrix size
nx=144;
ny=144;

% low saturation power and high saturation power
wl=0.5;
wh=1;

% RF frequency offset
x=[
    -10
    -8
    -6
    -5
    -4.75
    -4.5
    -4.25
    -4
    -3.75
    -3.5
    -3.25
    -3
    -2.75
    -2.5
    -2.25
    -2
    -1.75
    -1.5
    -1.25
    -1
    -0.75
    -0.5
    -0.25
    0
    0.25
    0.5
    0.75
    1
    1.25
    1.5
    1.75
    2
    2.25
    2.5
    2.75
    3
    3.25
    3.5
    3.75
    4
    4.25
    4.5
    4.75
    5
    6
    8
    10]*127.6596;


% RF frequency offset for LD fitting
    x_2pool=[
    -10
    -8
    -0.5
    -0.25
    0
    0.25
    0.5
    8
    10]*127.6596;

% RF freqency 
FreqArray = [-10,-8,-6,-5:0.25:5,6,8,10];
deta_cs =0.025;                                      
range_cs = -10:deta_cs:10;                             
interp1NOE = 801;                                  
referenceCriticalValue = 100;                     
cest1_0 = 0;         






% Sub1
sub_num=1;
% load ROI
load('Sub1\roi_brain')
load_nii_WM=load_nii('Sub1\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub1\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';

% T1 map
Imag1_T1=dicomread('Sub1\dcm_format\MR_1601_T1_MAP\9.dcm');

MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;

% read CEST data
%0p5 and 1p0
ls=dir('Sub1\dcm_format\MR_801_CEST_0_5uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_801_CEST_0_5uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p5uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p5uT(:,:,ii)=double(Imag1_0p5uT(:,:,ii))./mean(Imag1_0p5uT(:,:,[1,2,50,51]),3);
end


ls=dir('Sub1\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1uT(:,:,ii)=double(Imag1_1uT(:,:,ii))./mean(Imag1_1uT(:,:,[1,2,50,51]),3);
end

% CEST Z-spectrum without B0 correction
for ii=1:51
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
end

% no fat suppression was used in this paper. But this can be done by 
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p5uT_fatsuppression=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression=Zspectra_1uT;

% calculate DSP and then AREXDSP and MTRDSP without B0 correction
Zspectra_DPS_0p5uT_fatsuppression=1./(1+(1./Zspectra_1uT_fatsuppression-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression=1./Zspectra_0p5uT_fatsuppression-1./Zspectra_DPS_0p5uT_fatsuppression;
MTRDSP_0p5uT_fatsuppression=Zspectra_DPS_0p5uT_fatsuppression-Zspectra_0p5uT_fatsuppression;


% B0 correction 0.5uT
[Length,Width,Freq] = size(nImag1_0p5uT); 
nImag1_0p5uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_0p5uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_0p5uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

% B0 corretion 1uT
[Length,Width,Freq] = size(nImag1_1uT); 
nImag1_1uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_1uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_1uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

% residual MT corrected DSP maps with B0 correction
nImag1_DSP_0p5uT_corr=1./(1+(1./nImag1_1uT_corr-1).*((wl)/(wh))^2);
Image_MTR_APT_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,38)-nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;
Image_MTR_NOE_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,10)-nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;

Image_AREX_APT_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,38)-1./nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;
Image_AREX_NOE_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,10)-1./nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;


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

% R1
R1_map=1000./double(Imag1_T1);
R1_map(R1_map==inf)=nan;
R1_map(R1_map==-inf)=nan;
Value_R1_WM(sub_num)=nansum(nansum(R1_map.*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Value_R1_GM(sub_num)=nansum(nansum(R1_map.*roi_brain_GM))./nansum(nansum(roi_brain_GM));

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

Zspectra_0p5uT_WM_corr_single=Zspectra_0p5uT_WM_corr(:,sub_num);
% Lorentzian difference WM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound


Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_WM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_WM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_WM_corr(:,sub_num));
Zspectra_AREX_2pool_WM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_WM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_WM(sub_num);

Zspectra_0p5uT_GM_corr_single=Zspectra_0p5uT_GM_corr(:,sub_num);

% Lorentzian difference GM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_GM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_GM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_GM_corr(:,sub_num));
Zspectra_AREX_2pool_GM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_GM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_GM(sub_num);


Value_WM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_WM)))./nanmean(nanmean(roi_brain_WM));
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

Value_GM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_GM)))./nanmean(nanmean(roi_brain_GM));
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




% Sub2
sub_num=2;
load('Sub2\roi_brain')
load_nii_WM=load_nii('Sub2\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub2\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';
Imag1_T1=dicomread('Sub2\dcm_format\MR_1601_T1_MAP\9.dcm');

MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;


%0p5 and 1p0
ls=dir('Sub2\dcm_format\MR_801_CEST_0_5uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_801_CEST_0_5uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p5uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p5uT(:,:,ii)=double(Imag1_0p5uT(:,:,ii))./mean(Imag1_0p5uT(:,:,[1,2,50,51]),3);
end


ls=dir('Sub2\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1uT(:,:,ii)=double(Imag1_1uT(:,:,ii))./mean(Imag1_1uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
end

 Zspectra_0p5uT_fatsuppression=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression=1./(1+(1./Zspectra_1uT_fatsuppression-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression=1./Zspectra_0p5uT_fatsuppression-1./Zspectra_DPS_0p5uT_fatsuppression;
MTRDSP_0p5uT_fatsuppression=Zspectra_DPS_0p5uT_fatsuppression-Zspectra_0p5uT_fatsuppression;


% B0 correction 0.5uT
[Length,Width,Freq] = size(nImag1_0p5uT); 
nImag1_0p5uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_0p5uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_0p5uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

% B0 corretion 1uT
[Length,Width,Freq] = size(nImag1_1uT); 
nImag1_1uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_1uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_1uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

nImag1_DSP_0p5uT_corr=1./(1+(1./nImag1_1uT_corr-1).*((wl)/(wh))^2);
Image_MTR_APT_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,38)-nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;
Image_MTR_NOE_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,10)-nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;

Image_AREX_APT_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,38)-1./nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;
Image_AREX_NOE_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,10)-1./nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;

for ii=1:47
Zspectra_0p5uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_0p5uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
Zspectra_1uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_1uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end
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

% R1
R1_map=1000./double(Imag1_T1);
R1_map(R1_map==inf)=nan;
R1_map(R1_map==-inf)=nan;
Value_R1_WM(sub_num)=nansum(nansum(R1_map.*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Value_R1_GM(sub_num)=nansum(nansum(R1_map.*roi_brain_GM))./nansum(nansum(roi_brain_GM));

% asym
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

Zspectra_0p5uT_WM_corr_single=Zspectra_0p5uT_WM_corr(:,sub_num);
% Lorentzian difference WM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_WM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_WM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_WM_corr(:,sub_num));
Zspectra_AREX_2pool_WM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_WM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_WM(sub_num);

Zspectra_0p5uT_GM_corr_single=Zspectra_0p5uT_GM_corr(:,sub_num);
% Lorentzian difference GM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+38);
end


beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_GM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_GM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_GM_corr(:,sub_num));
Zspectra_AREX_2pool_GM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_GM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_GM(sub_num);



Value_WM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_WM)))./nanmean(nanmean(roi_brain_WM));
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

Value_GM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_GM)))./nanmean(nanmean(roi_brain_GM));
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





% Sub3
sub_num=3;
load('Sub3\roi_brain')
load_nii_WM=load_nii('Sub3\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub3\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';
Imag1_T1=dicomread('Sub3\dcm_format\MR_1601_T1_MAP\9.dcm');

MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;


%0p5 and 1p0
ls=dir('Sub3\dcm_format\MR_801_CEST_0_5uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_801_CEST_0_5uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p5uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p5uT(:,:,ii)=double(Imag1_0p5uT(:,:,ii))./mean(Imag1_0p5uT(:,:,[1,2,50,51]),3);
end


ls=dir('Sub3\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1uT(:,:,ii)=double(Imag1_1uT(:,:,ii))./mean(Imag1_1uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
end

 Zspectra_0p5uT_fatsuppression=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression=1./(1+(1./Zspectra_1uT_fatsuppression-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression=1./Zspectra_0p5uT_fatsuppression-1./Zspectra_DPS_0p5uT_fatsuppression;
MTRDSP_0p5uT_fatsuppression=Zspectra_DPS_0p5uT_fatsuppression-Zspectra_0p5uT_fatsuppression;


% B0 correction 0.5uT
[Length,Width,Freq] = size(nImag1_0p5uT); 
nImag1_0p5uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_0p5uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_0p5uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

% B0 corretion 1uT
[Length,Width,Freq] = size(nImag1_1uT); 
nImag1_1uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_1uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_1uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

nImag1_DSP_0p5uT_corr=1./(1+(1./nImag1_1uT_corr-1).*((wl)/(wh))^2);
Image_MTR_APT_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,38)-nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;
Image_MTR_NOE_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,10)-nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;

Image_AREX_APT_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,38)-1./nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;
Image_AREX_NOE_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,10)-1./nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;

for ii=1:47
Zspectra_0p5uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_0p5uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
Zspectra_1uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_1uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end
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

% R1
R1_map=1000./double(Imag1_T1);
R1_map(R1_map==inf)=nan;
R1_map(R1_map==-inf)=nan;
Value_R1_WM(sub_num)=nansum(nansum(R1_map.*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Value_R1_GM(sub_num)=nansum(nansum(R1_map.*roi_brain_GM))./nansum(nansum(roi_brain_GM));

% asym
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

Zspectra_0p5uT_WM_corr_single=Zspectra_0p5uT_WM_corr(:,sub_num);
% Lorentzian difference WM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_WM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_WM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_WM_corr(:,sub_num));
Zspectra_AREX_2pool_WM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_WM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_WM(sub_num);

Zspectra_0p5uT_GM_corr_single=Zspectra_0p5uT_GM_corr(:,sub_num);
% Lorentzian difference GM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

% Delta =[sep; R1S; R1W; R1M; R2M; mnotw];
Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_GM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_GM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_GM_corr(:,sub_num));
Zspectra_AREX_2pool_GM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_GM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_GM(sub_num);


Value_WM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_WM)))./nanmean(nanmean(roi_brain_WM));
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

Value_GM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_GM)))./nanmean(nanmean(roi_brain_GM));
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




% Sub4
sub_num=4;
load('Sub4\roi_brain')
load_nii_WM=load_nii('Sub4\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub4\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';
Imag1_T1=dicomread('Sub4\dcm_format\MR_1601_T1_MAP\9.dcm');

MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;


%0p5 and 1p0
ls=dir('Sub4\dcm_format\MR_801_CEST_0_5uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_801_CEST_0_5uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p5uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p5uT(:,:,ii)=double(Imag1_0p5uT(:,:,ii))./mean(Imag1_0p5uT(:,:,[1,2,50,51]),3);
end


ls=dir('Sub4\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1uT(:,:,ii)=double(Imag1_1uT(:,:,ii))./mean(Imag1_1uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
end

 Zspectra_0p5uT_fatsuppression=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression=1./(1+(1./Zspectra_1uT_fatsuppression-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression=1./Zspectra_0p5uT_fatsuppression-1./Zspectra_DPS_0p5uT_fatsuppression;
MTRDSP_0p5uT_fatsuppression=Zspectra_DPS_0p5uT_fatsuppression-Zspectra_0p5uT_fatsuppression;


% B0 correction 0.5uT
[Length,Width,Freq] = size(nImag1_0p5uT); 
nImag1_0p5uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_0p5uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_0p5uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

% B0 corretion 1uT
[Length,Width,Freq] = size(nImag1_1uT); 
nImag1_1uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_1uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_1uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

nImag1_DSP_0p5uT_corr=1./(1+(1./nImag1_1uT_corr-1).*((wl)/(wh))^2);
Image_MTR_APT_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,38)-nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;
Image_MTR_NOE_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,10)-nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;

Image_AREX_APT_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,38)-1./nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;
Image_AREX_NOE_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,10)-1./nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;

for ii=1:47
Zspectra_0p5uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_0p5uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
Zspectra_1uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_1uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end
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

% R1
R1_map=1000./double(Imag1_T1);
R1_map(R1_map==inf)=nan;
R1_map(R1_map==-inf)=nan;
Value_R1_WM(sub_num)=nansum(nansum(R1_map.*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Value_R1_GM(sub_num)=nansum(nansum(R1_map.*roi_brain_GM))./nansum(nansum(roi_brain_GM));

% asym
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

Zspectra_0p5uT_WM_corr_single=Zspectra_0p5uT_WM_corr(:,sub_num);
% Lorentzian difference WM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_WM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_WM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_WM_corr(:,sub_num));
Zspectra_AREX_2pool_WM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_WM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_WM(sub_num);

Zspectra_0p5uT_GM_corr_single=Zspectra_0p5uT_GM_corr(:,sub_num);
% Lorentzian difference GM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_GM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_GM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_GM_corr(:,sub_num));
Zspectra_AREX_2pool_GM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_GM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_GM(sub_num);



Value_WM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_WM)))./nanmean(nanmean(roi_brain_WM));
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

Value_GM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_GM)))./nanmean(nanmean(roi_brain_GM));
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




% Sub5
sub_num=5;
load('Sub5\roi_brain')
load_nii_WM=load_untouch_nii('Sub5\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_untouch_nii('Sub5\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';
Imag1_T1=dicomread('Sub5\dcm_format\MR_1601_T1_MAP\9.dcm');

MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;


%0p5 and 1p0
ls=dir('Sub5\dcm_format\MR_801_CEST_0_5uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_801_CEST_0_5uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p5uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p5uT(:,:,ii)=double(Imag1_0p5uT(:,:,ii))./mean(Imag1_0p5uT(:,:,[1,2,50,51]),3);
end


ls=dir('Sub5\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1uT(:,:,ii)=double(Imag1_1uT(:,:,ii))./mean(Imag1_1uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
end

 Zspectra_0p5uT_fatsuppression=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression=1./(1+(1./Zspectra_1uT_fatsuppression-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression=1./Zspectra_0p5uT_fatsuppression-1./Zspectra_DPS_0p5uT_fatsuppression;
MTRDSP_0p5uT_fatsuppression=Zspectra_DPS_0p5uT_fatsuppression-Zspectra_0p5uT_fatsuppression;


% B0 correction 0.5uT
[Length,Width,Freq] = size(nImag1_0p5uT); 
nImag1_0p5uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_0p5uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_0p5uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

% B0 corretion 1uT
[Length,Width,Freq] = size(nImag1_1uT); 
nImag1_1uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_1uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_1uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

nImag1_DSP_0p5uT_corr=1./(1+(1./nImag1_1uT_corr-1).*((wl)/(wh))^2);
Image_MTR_APT_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,38)-nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;
Image_MTR_NOE_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,10)-nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;

Image_AREX_APT_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,38)-1./nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;
Image_AREX_NOE_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,10)-1./nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;

for ii=1:47
Zspectra_0p5uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_0p5uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
Zspectra_1uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_1uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end
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

% R1
R1_map=1000./double(Imag1_T1);
R1_map(R1_map==inf)=nan;
R1_map(R1_map==-inf)=nan;
Value_R1_WM(sub_num)=nansum(nansum(R1_map.*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Value_R1_GM(sub_num)=nansum(nansum(R1_map.*roi_brain_GM))./nansum(nansum(roi_brain_GM));

% asym
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

Zspectra_0p5uT_WM_corr_single=Zspectra_0p5uT_WM_corr(:,sub_num);
% Lorentzian difference WM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_WM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_WM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_WM_corr(:,sub_num));
Zspectra_AREX_2pool_WM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_WM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_WM(sub_num);

Zspectra_0p5uT_GM_corr_single=Zspectra_0p5uT_GM_corr(:,sub_num);
% Lorentzian difference GM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_GM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_GM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_GM_corr(:,sub_num));
Zspectra_AREX_2pool_GM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_GM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_GM(sub_num);


Value_WM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_WM)))./nanmean(nanmean(roi_brain_WM));
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

Value_GM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_GM)))./nanmean(nanmean(roi_brain_GM));
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





% Sub6
sub_num=6;
load('Sub6\roi_brain')
load_nii_WM=load_nii('Sub6\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub6\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';
Imag1_T1=dicomread('Sub6\dcm_format\MR_1601_T1_MAP\9.dcm');

MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
MTR_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;
AREX_nImag1_0p5uT_corr_MT_corr(1:nx, 1:ny, 1:47)=0;


%0p5 and 1p0
ls=dir('Sub6\dcm_format\MR_801_CEST_0_5uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_801_CEST_0_5uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p5uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p5uT(:,:,ii)=double(Imag1_0p5uT(:,:,ii))./mean(Imag1_0p5uT(:,:,[1,2,50,51]),3);
end


ls=dir('Sub6\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_1201_CEST_1_0uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1uT(:,:,ii)=double(Imag1_1uT(:,:,ii))./mean(Imag1_1uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain))./nansum(nansum(roi_brain));
end

 Zspectra_0p5uT_fatsuppression=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression=1./(1+(1./Zspectra_1uT_fatsuppression-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression=1./Zspectra_0p5uT_fatsuppression-1./Zspectra_DPS_0p5uT_fatsuppression;
MTRDSP_0p5uT_fatsuppression=Zspectra_DPS_0p5uT_fatsuppression-Zspectra_0p5uT_fatsuppression;


% B0 correction 0.5uT
[Length,Width,Freq] = size(nImag1_0p5uT); 
nImag1_0p5uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_0p5uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_0p5uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

% B0 corretion 1uT
[Length,Width,Freq] = size(nImag1_1uT); 
nImag1_1uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_1uT(xx,yy,3:49).*roi_brain(xx,yy))';
            if Raw_Z(1)==0 | isnan(Raw_Z(1))==1
                  Mask1(xx,yy,:)=0;
            else
                  Mask1(xx,yy,:)=1;
                AA = interp1(FreqArray,Raw_Z,range_cs,'spline');
                InterplCESTData(xx,yy,:) = AA;
            end
    end
end
    

WASSRmap = zeros(Length,Width);                      
Corrected_Z = zeros(Length,Width,interp1NOE);          
WASSRAcquisitionNumber = interp1NOE;
detacenter = 0;

for xx = 1:Length
    for yy = 1:Width
            WASSRvalue = squeeze(InterplCESTData(xx,yy,:));
            [minvalue,position] = min(WASSRvalue);
            if position > 1 && position < WASSRAcquisitionNumber
                deta = WASSRvalue(position-1) - WASSRvalue(position+1);
                if deta > 0
                    detacenter = 0.5*deta_cs -  (0.5*deta_cs)/((WASSRvalue(position-1) - minvalue)/(WASSRvalue(position+1) - minvalue));
                elseif deta < 0
                    detacenter = -0.5*deta_cs +  (0.5*deta_cs)/((WASSRvalue(position+1) - minvalue)/(WASSRvalue(position-1) - minvalue));
                elseif deta == 0
                    detacenter = 0;
                end
            end
            center = ((WASSRAcquisitionNumber+1) / 2 - position)*deta_cs-detacenter;
            WASSRmap(xx,yy) = center;
            detacenter = 0;
            
            corrected_cs =range_cs + double(center);        
            Uncorrected_Z = squeeze(InterplCESTData(xx,yy,:));
            if Uncorrected_Z(1)>0.05 && abs(center) < 2     
                Corrected_Z(xx,yy,:) = interp1(corrected_cs,Uncorrected_Z, range_cs,'spline');
                TarValue_2 = interp1(corrected_cs,Uncorrected_Z, FreqArray,'spline');
                nImag1_1uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

nImag1_DSP_0p5uT_corr=1./(1+(1./nImag1_1uT_corr-1).*((wl)/(wh))^2);
Image_MTR_APT_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,38)-nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;
Image_MTR_NOE_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,10)-nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_brain;

Image_AREX_APT_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,38)-1./nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;
Image_AREX_NOE_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,10)-1./nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_brain;

for ii=1:47
Zspectra_0p5uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_0p5uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
Zspectra_1uT_WM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Zspectra_1uT_GM_corr(ii,sub_num)=nansum(nansum(nImag1_1uT_corr(:,:,ii).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end
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

% R1
R1_map=1000./double(Imag1_T1);
R1_map(R1_map==inf)=nan;
R1_map(R1_map==-inf)=nan;
Value_R1_WM(sub_num)=nansum(nansum(R1_map.*roi_brain_WM))./nansum(nansum(roi_brain_WM));
Value_R1_GM(sub_num)=nansum(nansum(R1_map.*roi_brain_GM))./nansum(nansum(roi_brain_GM));

% asym
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

Zspectra_0p5uT_WM_corr_single=Zspectra_0p5uT_WM_corr(:,sub_num);
% Lorentzian difference WM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_WM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_WM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_WM_corr(:,sub_num));
Zspectra_AREX_2pool_WM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_WM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_WM(sub_num);

Zspectra_0p5uT_WM_corr_single=Zspectra_0p5uT_WM_corr(:,sub_num);
% Lorentzian difference WM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_WM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_WM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_WM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_WM_corr(:,sub_num));
Zspectra_AREX_2pool_WM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_WM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_WM(sub_num);

Zspectra_0p5uT_GM_corr_single=Zspectra_0p5uT_GM_corr(:,sub_num);
% Lorentzian difference GM
for i=1:1:2
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i);
end
for i=3:1:7
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+19);
end
for i=8:1:9
    sig_2pool(i)=1-Zspectra_0p5uT_GM_corr_single(i+38);
end

beta0_2pool= [0.9,       0,       1.4*127.66,      0.1,   0*127.66         25*127.66]; % initial test
lb_2pool=    [0.02,  -1*127.66,   0.1*127.66,      0,    -4*127.66,     10*127.66]; % lower bound
ub_2pool=    [1,     1*127.66,    10*127.66,       1,     4*127.66,     100*127.66]; % upper bound

Delta=[1]; % constants

options=optimset('lsqcurvefit') ; 
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

[beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool', sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

sig_simu_2pool=matsolv_2pool(beta_2pool,x,Delta);

Zspectra_reference_2pool_GM(:,sub_num)=1-sig_simu_2pool;
Zspectra_MTR_2pool_GM(:,sub_num)=1-sig_simu_2pool-squeeze(Zspectra_0p5uT_GM_corr(:,sub_num));
Zspectra_AREX_2pool_GM(:,sub_num)=(1./squeeze(Zspectra_0p5uT_GM_corr(:,sub_num))-1./(1-sig_simu_2pool))*Value_R1_GM(sub_num);

Value_WM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_WM)))./nanmean(nanmean(roi_brain_WM));
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

Value_GM_anatomy(sub_num)=nanmean(nanmean((nanmean(Imag1_0p5uT(:,:,[1,2,50,51]),3).*roi_brain_GM)))./nanmean(nanmean(roi_brain_GM));
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


sub_range=[1:6]
% Fig. 6a
figure (1)
hold on
plot(flip(nanmean(Zspectra_0p5uT_WM_corr(:,sub_range),2)));
plot(flip(nanmean(Zspectra_1uT_WM_corr(:,sub_range),2)));

% Fig. 6c
figure (2)
hold on
plot(flip(nanmean(Zspectra_0p5uT_WM_corr(:,sub_range),2)));
plot(flip(nanmean(Zspectra_DSP_0p5uT_WM_corr(:,sub_range),2)));

% Fig. 6e
figure (3)
hold on
plot(flip(nanmean(Zspectra_MTR_WM_corr(:,sub_range)')));
%errorbar(flip(nanmean(Zspectra_MTR_WM_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_WM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_MTR_WM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_WM_corr_MT_corr(:,sub_range)')));

% Fig. 6g
figure (4)
hold on
plot(flip(nanmean(Zspectra_AREX_WM_corr(:,sub_range)')));
%errorbar(flip(nanmean(Zspectra_AREX_WM_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_WM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_AREX_WM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_WM_corr_MT_corr(:,sub_range)')));


% Fig. 6b
figure (1)
hold on
plot(flip(nanmean(Zspectra_0p5uT_GM_corr(:,sub_range),2)));
plot(flip(nanmean(Zspectra_1uT_GM_corr(:,sub_range),2)));

% Fig. 6d
figure (2)
hold on
plot(flip(nanmean(Zspectra_0p5uT_GM_corr(:,sub_range),2)));
plot(flip(nanmean(Zspectra_DSP_0p5uT_GM_corr(:,sub_range),2)));

% Fig. 6f
figure (3)
hold on
plot(flip(nanmean(Zspectra_MTR_GM_corr(:,sub_range)')));
%errorbar(flip(nanmean(Zspectra_MTR_GM_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_GM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_MTR_GM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_GM_corr_MT_corr(:,sub_range)')));

% Fig. 6h
figure (4)
hold on
plot(flip(nanmean(Zspectra_AREX_GM_corr(:,sub_range)')));
%errorbar(flip(nanmean(Zspectra_AREX_GM_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_GM_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_AREX_GM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_GM_corr_MT_corr(:,sub_range)')));









% Fig. 7a
figure (1)
hold on
plot(flip(nanmean(Zspectra_0p5uT_WM_corr(:,sub_range),2)));
plot(flip(nanmean(Zspectra_reference_2pool_WM(:,sub_range),2)));

% Fig. 7c
figure (2)
hold on
errorbar(flip(nanmean(Zspectra_MTR_2pool_WM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_WM(:,sub_range)')));

% Fig. 7e
figure (3)
hold on
errorbar(flip(nanmean(Zspectra_AREX_2pool_WM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_WM(:,sub_range)')));

% Fig. 7g
figure (4)
hold on
errorbar((nanmean(fMTRasym_0p5uT_WM_corr(:,sub_range)')), (nanstd(fMTRasym_0p5uT_WM_corr(:,sub_range)')));

% Fig. 7i
figure (5)
hold on
errorbar((nanmean(fAREXasym_0p5uT_WM_corr(:,sub_range)')), (nanstd(fAREXasym_0p5uT_WM_corr(:,sub_range)')));



% Fig. 7b
figure (1)
hold on
plot(flip(nanmean(Zspectra_0p5uT_GM_corr(:,sub_range),2)));
plot(flip(nanmean(Zspectra_reference_2pool_GM(:,sub_range),2)));

% Fig. 7d
figure (2)
hold on
errorbar(flip(nanmean(Zspectra_MTR_2pool_GM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_GM(:,sub_range)')));

% Fig. 7f
figure (3)
hold on
errorbar(flip(nanmean(Zspectra_AREX_2pool_GM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_GM(:,sub_range)')));

% Fig. 7h
figure (4)
hold on
errorbar((nanmean(fMTRasym_0p5uT_GM_corr(:,sub_range)')), (nanstd(fMTRasym_0p5uT_GM_corr(:,sub_range)')));

% Fig. 7j
figure (5)
hold on
errorbar((nanmean(fAREXasym_0p5uT_GM_corr(:,sub_range)')), (nanstd(fAREXasym_0p5uT_GM_corr(:,sub_range)')));



%Fig.8a
figure (1)
hold on
errorbar(flip(nanmean(Zspectra_MTR_WM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_WM_corr_MT_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_MTR_2pool_WM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_WM(:,sub_range)')));
errorbar((nanmean(fMTRasym_0p5uT_WM_corr(:,sub_range)')), (nanstd(fMTRasym_0p5uT_WM_corr(:,sub_range)')));

%Fig.8b
figure (2)
hold on
errorbar(flip(nanmean(Zspectra_MTR_GM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_MTR_GM_corr_MT_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_MTR_2pool_GM(:,sub_range)')), flip(nanstd(Zspectra_MTR_2pool_GM(:,sub_range)')));
errorbar((nanmean(fMTRasym_0p5uT_GM_corr(:,sub_range)')), (nanstd(fMTRasym_0p5uT_GM_corr(:,sub_range)')));

%Fig.8c
figure (3)
hold on
errorbar(flip(nanmean(Zspectra_AREX_WM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_WM_corr_MT_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_AREX_2pool_WM(:,sub_range)')), flip(nanstd(Zspectra_AREX_2pool_WM(:,sub_range)')));
errorbar((nanmean(fAREXasym_0p5uT_WM_corr(:,sub_range)')), (nanstd(fAREXasym_0p5uT_WM_corr(:,sub_range)')));

%Fig.8d
figure (4)
hold on
errorbar(flip(nanmean(Zspectra_AREX_GM_corr_MT_corr(:,sub_range)')), flip(nanstd(Zspectra_AREX_GM_corr_MT_corr(:,sub_range)')));
errorbar(flip(nanmean(Zspectra_AREX_2pool_GM(:,sub_range)')), flip(nanstd(Zspectra_AREX_2pool_GM(:,sub_range)')));
errorbar((nanmean(fAREXasym_0p5uT_GM_corr(:,sub_range)')), (nanstd(fAREXasym_0p5uT_GM_corr(:,sub_range)')));







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

figure(1);
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

figure(1);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_GM_anatomy(sub_range), Value_WM_anatomy(sub_range))
ylabel('Anatomy (a.u.)')
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

figure(2);
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

figure(2);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_R1_WM(sub_range), Value_R1_GM(sub_range))
ylabel('R_{1obs} (s^{-1})')
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

figure(3);
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

figure(3);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_GM_MTR_corr_APT(sub_range), Value_WM_MTR_corr_APT(sub_range))
ylabel('MTR_{DSP} at 3.5ppm (%)')
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

figure(4);
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

figure(4);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_AREX_corr_APT(sub_range), Value_GM_AREX_corr_APT(sub_range))
ylabel('AREX_{DSP} at 3.5ppm (%s^{-1})')
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

figure(5);
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

figure(5);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_MTR_2pool_APT(sub_range), Value_GM_MTR_2pool_APT(sub_range))
ylabel('MTR_{LD} at 3.5ppm (%)')
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

figure(6);
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

figure(6);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_AREX_2pool_APT(sub_range), Value_GM_AREX_2pool_APT(sub_range))
ylabel('AREX_{LD} at 3.5ppm (%s^{-1})')
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

figure(7);
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

figure(7);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_MTRasym_corr(sub_range), Value_GM_MTRasym_corr(sub_range))
ylabel('MTR_{asym} at 3.5ppm (%)')
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

figure(8);
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

figure(8);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_AREXasym_corr(sub_range), Value_GM_AREXasym_corr(sub_range))
ylabel('AREX_{asym} at 3.5ppm (%s^{-1})')
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

figure(9);
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

figure(9);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_MTR_corr_NOE(sub_range), Value_GM_MTR_corr_NOE(sub_range))
ylabel('MTR_{DSP} at -3.5ppm (%)')
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

figure(10);
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

figure(10);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_AREX_corr_NOE(sub_range), Value_GM_AREX_corr_NOE(sub_range))
ylabel('AREX_{DSP} at -3.5ppm (%s^{-1})')
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

figure(11);
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

figure(11);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_MTR_2pool_NOE(sub_range), Value_GM_MTR_2pool_NOE(sub_range))
ylabel('MTR_{LD} at -3.5ppm (%)')
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

figure(12);
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

figure(12);
hold on;
bar(x_zhu,y_zhu);
errorbar(x_zhu,y_zhu,L_zhu,U_zhu,'Marker','none','LineStyle','none')

[a,b]=ttest(Value_WM_AREX_2pool_APT(sub_range), Value_GM_AREX_2pool_APT(sub_range))
ylabel('AREX_{LD} at -3.5ppm (%s^{-1})')
path = sprintf('%s%d.fig','FigS7\',12);
set(gca,'xtick',[]);set(gca,'ytick',[]);




