clear all; close all; clc;

% matrix size
nx=144;
ny=144;

% low saturation power and high saturation power
wl=0.5;
wh=1;


% RF frequency offset
xx=[100
    -100
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
    10
    100
    -100];




% optimizing WM
% Sub1
rat_num=1;
% load ROI
load('Sub1\roi_brain')
load_nii_WM=load_nii('Sub1\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub1\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';

% read CEST data
%0p3uT and 0p6uT
ls=dir('Sub1\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub1\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression=1./(1+(1./Zspectra_0p6uT_fatsuppression-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression(:,rat_num)=1./Zspectra_0p3uT_fatsuppression-1./Zspectra_DPS_0p3uT_fatsuppression;
MTRDSP_0p3uT_fatsuppression(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression-Zspectra_0p3uT_fatsuppression;







%0p4uT and 0p8uT
ls=dir('Sub1\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub1\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_WM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_WM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_WM=1./(1+(1./Zspectra_0p8uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_WM-1./Zspectra_DPS_0p4uT_fatsuppression_WM;
MTRDSP_0p4uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_WM-Zspectra_0p4uT_fatsuppression_WM;









%0p5uT and 1p0uT
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

for ii=1:51
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_WM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_WM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_WM=1./(1+(1./Zspectra_1uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_WM-1./Zspectra_DPS_0p5uT_fatsuppression_WM;
MTRDSP_0p5uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_WM-Zspectra_0p5uT_fatsuppression_WM;






%0p6uT and 1p2uT
ls=dir('Sub1\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub1\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_WM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_WM=1./(1+(1./Zspectra_1p2uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_WM-1./Zspectra_DPS_0p6uT_fatsuppression_WM;
MTRDSP_0p6uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_WM-Zspectra_0p6uT_fatsuppression_WM;





%0p7 and 1p4
ls=dir('Sub1\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub1\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression_WM=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression_WM=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_WM=1./(1+(1./Zspectra_1p4uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_WM-1./Zspectra_DPS_0p7uT_fatsuppression_WM;
MTRDSP_0p7uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_WM-Zspectra_0p7uT_fatsuppression_WM;











% Sub2
rat_num=2;
load('Sub2\roi_brain')
load_nii_WM=load_nii('Sub2\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub2\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';


%0p3uT and 0p6uT
ls=dir('Sub2\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub2\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression_WM=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression_WM=1./(1+(1./Zspectra_0p6uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p3uT_fatsuppression_WM-1./Zspectra_DPS_0p3uT_fatsuppression_WM;
MTRDSP_0p3uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression_WM-Zspectra_0p3uT_fatsuppression_WM;







%0p4uT and 0p8uT
ls=dir('Sub2\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub2\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_WM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_WM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_WM=1./(1+(1./Zspectra_0p8uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_WM-1./Zspectra_DPS_0p4uT_fatsuppression_WM;
MTRDSP_0p4uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_WM-Zspectra_0p4uT_fatsuppression_WM;









%0p5uT and 1p0uT
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
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_WM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_WM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_WM=1./(1+(1./Zspectra_1uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_WM-1./Zspectra_DPS_0p5uT_fatsuppression_WM;
MTRDSP_0p5uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_WM-Zspectra_0p5uT_fatsuppression_WM;






%0p6uT and 1p2uT
ls=dir('Sub2\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub2\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_WM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_WM=1./(1+(1./Zspectra_1p2uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_WM-1./Zspectra_DPS_0p6uT_fatsuppression_WM;
MTRDSP_0p6uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_WM-Zspectra_0p6uT_fatsuppression_WM;





%0p7 and 1p4
ls=dir('Sub2\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub2\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression_WM=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression_WM=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_WM=1./(1+(1./Zspectra_1p4uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_WM-1./Zspectra_DPS_0p7uT_fatsuppression_WM;
MTRDSP_0p7uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_WM-Zspectra_0p7uT_fatsuppression_WM;










% Sub3
rat_num=3;
load('Sub3\roi_brain')
load_nii_WM=load_nii('Sub3\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub3\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';


%0p3uT and 0p6uT
ls=dir('Sub3\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub3\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression_WM=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression_WM=1./(1+(1./Zspectra_0p6uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p3uT_fatsuppression_WM-1./Zspectra_DPS_0p3uT_fatsuppression_WM;
MTRDSP_0p3uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression_WM-Zspectra_0p3uT_fatsuppression_WM;







%0p4uT and 0p8uT
ls=dir('Sub3\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub3\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_WM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_WM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_WM=1./(1+(1./Zspectra_0p8uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_WM-1./Zspectra_DPS_0p4uT_fatsuppression_WM;
MTRDSP_0p4uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_WM-Zspectra_0p4uT_fatsuppression_WM;









%0p5uT and 1p0uT
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
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_WM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_WM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_WM=1./(1+(1./Zspectra_1uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_WM-1./Zspectra_DPS_0p5uT_fatsuppression_WM;
MTRDSP_0p5uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_WM-Zspectra_0p5uT_fatsuppression_WM;






%0p6uT and 1p2uT
ls=dir('Sub3\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub3\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_WM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_WM=1./(1+(1./Zspectra_1p2uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_WM-1./Zspectra_DPS_0p6uT_fatsuppression_WM;
MTRDSP_0p6uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_WM-Zspectra_0p6uT_fatsuppression_WM;





%0p7 and 1p4
ls=dir('Sub3\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub3\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_WM=1./(1+(1./Zspectra_1p4uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_WM-1./Zspectra_DPS_0p7uT_fatsuppression_WM;
MTRDSP_0p7uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_WM-Zspectra_0p7uT_fatsuppression_WM;








% Sub4
rat_num=4;
load('Sub4\roi_brain')
load_nii_WM=load_nii('Sub4\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub4\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';


%0p3uT and 0p6uT
ls=dir('Sub4\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub4\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression_WM=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression_WM=1./(1+(1./Zspectra_0p6uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p3uT_fatsuppression_WM-1./Zspectra_DPS_0p3uT_fatsuppression_WM;
MTRDSP_0p3uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression_WM-Zspectra_0p3uT_fatsuppression_WM;







%0p4uT and 0p8uT
ls=dir('Sub4\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub4\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_WM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_WM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_WM=1./(1+(1./Zspectra_0p8uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_WM-1./Zspectra_DPS_0p4uT_fatsuppression_WM;
MTRDSP_0p4uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_WM-Zspectra_0p4uT_fatsuppression_WM;









%0p5uT and 1p0uT
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
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_WM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_WM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_WM=1./(1+(1./Zspectra_1uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_WM-1./Zspectra_DPS_0p5uT_fatsuppression_WM;
MTRDSP_0p5uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_WM-Zspectra_0p5uT_fatsuppression_WM;






%0p6uT and 1p2uT
ls=dir('Sub4\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub4\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_WM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_WM=1./(1+(1./Zspectra_1p2uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_WM-1./Zspectra_DPS_0p6uT_fatsuppression_WM;
MTRDSP_0p6uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_WM-Zspectra_0p6uT_fatsuppression_WM;





%0p7 and 1p4
ls=dir('Sub4\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub4\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression_WM=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression_WM=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_WM=1./(1+(1./Zspectra_1p4uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_WM-1./Zspectra_DPS_0p7uT_fatsuppression_WM;
MTRDSP_0p7uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_WM-Zspectra_0p7uT_fatsuppression_WM;








% Sub5
rat_num=5;
load('Sub5\roi_brain')
load_nii_WM=load_untouch_nii('Sub5\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_untouch_nii('Sub5\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';


%0p3uT and 0p6uT
ls=dir('Sub5\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub5\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression_WM=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression_WM=1./(1+(1./Zspectra_0p6uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p3uT_fatsuppression_WM-1./Zspectra_DPS_0p3uT_fatsuppression_WM;
MTRDSP_0p3uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression_WM-Zspectra_0p3uT_fatsuppression_WM;







%0p4uT and 0p8uT
ls=dir('Sub5\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub5\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_WM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_WM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_WM=1./(1+(1./Zspectra_0p8uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_WM-1./Zspectra_DPS_0p4uT_fatsuppression_WM;
MTRDSP_0p4uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_WM-Zspectra_0p4uT_fatsuppression_WM;









%0p5uT and 1p0uT
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
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_WM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_WM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_WM=1./(1+(1./Zspectra_1uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_WM-1./Zspectra_DPS_0p5uT_fatsuppression_WM;
MTRDSP_0p5uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_WM-Zspectra_0p5uT_fatsuppression_WM;






%0p6uT and 1p2uT
ls=dir('Sub5\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub5\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_WM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_WM=1./(1+(1./Zspectra_1p2uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_WM-1./Zspectra_DPS_0p6uT_fatsuppression_WM;
MTRDSP_0p6uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_WM-Zspectra_0p6uT_fatsuppression_WM;





%0p7 and 1p4
ls=dir('Sub5\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub5\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression_WM=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression_WM=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_WM=1./(1+(1./Zspectra_1p4uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_WM-1./Zspectra_DPS_0p7uT_fatsuppression_WM;
MTRDSP_0p7uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_WM-Zspectra_0p7uT_fatsuppression_WM;









% Sub6
rat_num=6;
load('Sub6\roi_brain')
load_nii_WM=load_nii('Sub6\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub6\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';

%0p3uT and 0p6uT
ls=dir('Sub6\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub6\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression_WM=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression_WM=1./(1+(1./Zspectra_0p6uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p3uT_fatsuppression_WM-1./Zspectra_DPS_0p3uT_fatsuppression_WM;
MTRDSP_0p3uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression_WM-Zspectra_0p3uT_fatsuppression_WM;







%0p4uT and 0p8uT
ls=dir('Sub6\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub6\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_WM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_WM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_WM=1./(1+(1./Zspectra_0p8uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_WM-1./Zspectra_DPS_0p4uT_fatsuppression_WM;
MTRDSP_0p4uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_WM-Zspectra_0p4uT_fatsuppression_WM;









%0p5uT and 1p0uT
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
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_WM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_WM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_WM=1./(1+(1./Zspectra_1uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_WM-1./Zspectra_DPS_0p5uT_fatsuppression_WM;
MTRDSP_0p5uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_WM-Zspectra_0p5uT_fatsuppression_WM;






%0p6uT and 1p2uT
ls=dir('Sub6\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub6\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_WM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_WM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_WM=1./(1+(1./Zspectra_1p2uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_WM-1./Zspectra_DPS_0p6uT_fatsuppression_WM;
MTRDSP_0p6uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_WM-Zspectra_0p6uT_fatsuppression_WM;





%0p7 and 1p4
ls=dir('Sub6\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub6\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_WM))./nansum(nansum(roi_brain_WM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression_WM=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression_WM=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_WM=1./(1+(1./Zspectra_1p4uT_fatsuppression_WM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_WM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_WM-1./Zspectra_DPS_0p7uT_fatsuppression_WM;
MTRDSP_0p7uT_fatsuppression_WM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_WM-Zspectra_0p7uT_fatsuppression_WM;







% optimizing GM
% Sub1
rat_num=1;
% load ROI
load('Sub1\roi_brain')
load_nii_WM=load_nii('Sub1\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub1\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';

% read CEST data
%0p3uT and 0p6uT
ls=dir('Sub1\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub1\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression=1./(1+(1./Zspectra_0p6uT_fatsuppression-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression(:,rat_num)=1./Zspectra_0p3uT_fatsuppression-1./Zspectra_DPS_0p3uT_fatsuppression;
MTRDSP_0p3uT_fatsuppression(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression-Zspectra_0p3uT_fatsuppression;







%0p4uT and 0p8uT
ls=dir('Sub1\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub1\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_GM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_GM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_GM=1./(1+(1./Zspectra_0p8uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_GM-1./Zspectra_DPS_0p4uT_fatsuppression_GM;
MTRDSP_0p4uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_GM-Zspectra_0p4uT_fatsuppression_GM;









%0p5uT and 1p0uT
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

for ii=1:51
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_GM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_GM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_GM=1./(1+(1./Zspectra_1uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_GM-1./Zspectra_DPS_0p5uT_fatsuppression_GM;
MTRDSP_0p5uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_GM-Zspectra_0p5uT_fatsuppression_GM;






%0p6uT and 1p2uT
ls=dir('Sub1\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub1\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_GM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_GM=1./(1+(1./Zspectra_1p2uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_GM-1./Zspectra_DPS_0p6uT_fatsuppression_GM;
MTRDSP_0p6uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_GM-Zspectra_0p6uT_fatsuppression_GM;





%0p7 and 1p4
ls=dir('Sub1\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub1\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub1\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression_GM=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression_GM=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_GM=1./(1+(1./Zspectra_1p4uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_GM-1./Zspectra_DPS_0p7uT_fatsuppression_GM;
MTRDSP_0p7uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_GM-Zspectra_0p7uT_fatsuppression_GM;











% Sub2
rat_num=2;
load('Sub2\roi_brain')
load_nii_WM=load_nii('Sub2\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub2\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';


%0p3uT and 0p6uT
ls=dir('Sub2\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub2\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression_GM=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression_GM=1./(1+(1./Zspectra_0p6uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p3uT_fatsuppression_GM-1./Zspectra_DPS_0p3uT_fatsuppression_GM;
MTRDSP_0p3uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression_GM-Zspectra_0p3uT_fatsuppression_GM;







%0p4uT and 0p8uT
ls=dir('Sub2\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub2\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_GM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_GM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_GM=1./(1+(1./Zspectra_0p8uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_GM-1./Zspectra_DPS_0p4uT_fatsuppression_GM;
MTRDSP_0p4uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_GM-Zspectra_0p4uT_fatsuppression_GM;









%0p5uT and 1p0uT
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
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_GM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_GM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_GM=1./(1+(1./Zspectra_1uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_GM-1./Zspectra_DPS_0p5uT_fatsuppression_GM;
MTRDSP_0p5uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_GM-Zspectra_0p5uT_fatsuppression_GM;






%0p6uT and 1p2uT
ls=dir('Sub2\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub2\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_GM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_GM=1./(1+(1./Zspectra_1p2uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_GM-1./Zspectra_DPS_0p6uT_fatsuppression_GM;
MTRDSP_0p6uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_GM-Zspectra_0p6uT_fatsuppression_GM;





%0p7 and 1p4
ls=dir('Sub2\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub2\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub2\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression_GM=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression_GM=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_GM=1./(1+(1./Zspectra_1p4uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_GM-1./Zspectra_DPS_0p7uT_fatsuppression_GM;
MTRDSP_0p7uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_GM-Zspectra_0p7uT_fatsuppression_GM;










% Sub3
rat_num=3;
load('Sub3\roi_brain')
load_nii_WM=load_nii('Sub3\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub3\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';


%0p3uT and 0p6uT
ls=dir('Sub3\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub3\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression_GM=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression_GM=1./(1+(1./Zspectra_0p6uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p3uT_fatsuppression_GM-1./Zspectra_DPS_0p3uT_fatsuppression_GM;
MTRDSP_0p3uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression_GM-Zspectra_0p3uT_fatsuppression_GM;







%0p4uT and 0p8uT
ls=dir('Sub3\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub3\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_GM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_GM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_GM=1./(1+(1./Zspectra_0p8uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_GM-1./Zspectra_DPS_0p4uT_fatsuppression_GM;
MTRDSP_0p4uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_GM-Zspectra_0p4uT_fatsuppression_GM;









%0p5uT and 1p0uT
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
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_GM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_GM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_GM=1./(1+(1./Zspectra_1uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_GM-1./Zspectra_DPS_0p5uT_fatsuppression_GM;
MTRDSP_0p5uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_GM-Zspectra_0p5uT_fatsuppression_GM;






%0p6uT and 1p2uT
ls=dir('Sub3\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub3\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_GM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_GM=1./(1+(1./Zspectra_1p2uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_GM-1./Zspectra_DPS_0p6uT_fatsuppression_GM;
MTRDSP_0p6uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_GM-Zspectra_0p6uT_fatsuppression_GM;





%0p7 and 1p4
ls=dir('Sub3\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub3\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub3\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_GM=1./(1+(1./Zspectra_1p4uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_GM-1./Zspectra_DPS_0p7uT_fatsuppression_GM;
MTRDSP_0p7uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_GM-Zspectra_0p7uT_fatsuppression_GM;








% Sub4
rat_num=4;
load('Sub4\roi_brain')
load_nii_WM=load_nii('Sub4\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub4\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';


%0p3uT and 0p6uT
ls=dir('Sub4\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub4\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression_GM=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression_GM=1./(1+(1./Zspectra_0p6uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p3uT_fatsuppression_GM-1./Zspectra_DPS_0p3uT_fatsuppression_GM;
MTRDSP_0p3uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression_GM-Zspectra_0p3uT_fatsuppression_GM;







%0p4uT and 0p8uT
ls=dir('Sub4\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub4\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_GM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_GM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_GM=1./(1+(1./Zspectra_0p8uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_GM-1./Zspectra_DPS_0p4uT_fatsuppression_GM;
MTRDSP_0p4uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_GM-Zspectra_0p4uT_fatsuppression_GM;









%0p5uT and 1p0uT
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
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_GM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_GM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_GM=1./(1+(1./Zspectra_1uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_GM-1./Zspectra_DPS_0p5uT_fatsuppression_GM;
MTRDSP_0p5uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_GM-Zspectra_0p5uT_fatsuppression_GM;






%0p6uT and 1p2uT
ls=dir('Sub4\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub4\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_GM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_GM=1./(1+(1./Zspectra_1p2uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_GM-1./Zspectra_DPS_0p6uT_fatsuppression_GM;
MTRDSP_0p6uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_GM-Zspectra_0p6uT_fatsuppression_GM;





%0p7 and 1p4
ls=dir('Sub4\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub4\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub4\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression_GM=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression_GM=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_GM=1./(1+(1./Zspectra_1p4uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_GM-1./Zspectra_DPS_0p7uT_fatsuppression_GM;
MTRDSP_0p7uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_GM-Zspectra_0p7uT_fatsuppression_GM;








% Sub5
rat_num=5;
load('Sub5\roi_brain')
load_nii_WM=load_untouch_nii('Sub5\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_untouch_nii('Sub5\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';


%0p3uT and 0p6uT
ls=dir('Sub5\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub5\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression_GM=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression_GM=1./(1+(1./Zspectra_0p6uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p3uT_fatsuppression_GM-1./Zspectra_DPS_0p3uT_fatsuppression_GM;
MTRDSP_0p3uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression_GM-Zspectra_0p3uT_fatsuppression_GM;







%0p4uT and 0p8uT
ls=dir('Sub5\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub5\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_GM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_GM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_GM=1./(1+(1./Zspectra_0p8uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_GM-1./Zspectra_DPS_0p4uT_fatsuppression_GM;
MTRDSP_0p4uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_GM-Zspectra_0p4uT_fatsuppression_GM;









%0p5uT and 1p0uT
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
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_GM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_GM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_GM=1./(1+(1./Zspectra_1uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_GM-1./Zspectra_DPS_0p5uT_fatsuppression_GM;
MTRDSP_0p5uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_GM-Zspectra_0p5uT_fatsuppression_GM;






%0p6uT and 1p2uT
ls=dir('Sub5\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub5\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_GM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_GM=1./(1+(1./Zspectra_1p2uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_GM-1./Zspectra_DPS_0p6uT_fatsuppression_GM;
MTRDSP_0p6uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_GM-Zspectra_0p6uT_fatsuppression_GM;





%0p7 and 1p4
ls=dir('Sub5\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub5\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub5\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression_GM=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression_GM=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_GM=1./(1+(1./Zspectra_1p4uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_GM-1./Zspectra_DPS_0p7uT_fatsuppression_GM;
MTRDSP_0p7uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_GM-Zspectra_0p7uT_fatsuppression_GM;









% Sub6
rat_num=6;
load('Sub6\roi_brain')
load_nii_WM=load_nii('Sub6\dcm_format\WM_mask.nii');
roi_brain_WM= fliplr(flipud(double(load_nii_WM.img)))';
load_nii_GM=load_nii('Sub6\dcm_format\GM_mask.nii');
roi_brain_GM= fliplr(flipud(double(load_nii_GM.img)))';

%0p3uT and 0p6uT
ls=dir('Sub6\dcm_format\MR_601_CEST_0_3uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_601_CEST_0_3uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p3uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p3uT(:,:,ii)=double(Imag1_0p3uT(:,:,ii))./mean(Imag1_0p3uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub6\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p3uT(ii)=nansum(nansum(double(nImag1_0p3uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p3uT_fatsuppression=(Zspectra_0p3uT(:)-Zspectra_0p3uT(26))./(1-Zspectra_0p3uT(26))
%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
 Zspectra_0p3uT_fatsuppression_GM=Zspectra_0p3uT;
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;

Zspectra_DPS_0p3uT_fatsuppression_GM=1./(1+(1./Zspectra_0p6uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p3uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p3uT_fatsuppression_GM-1./Zspectra_DPS_0p3uT_fatsuppression_GM;
MTRDSP_0p3uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p3uT_fatsuppression_GM-Zspectra_0p3uT_fatsuppression_GM;







%0p4uT and 0p8uT
ls=dir('Sub6\dcm_format\MR_701_CEST_0_4uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_701_CEST_0_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p4uT(:,:,ii)=double(Imag1_0p4uT(:,:,ii))./mean(Imag1_0p4uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub6\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_1101_CEST_0_8uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p8uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p8uT(:,:,ii)=double(Imag1_0p8uT(:,:,ii))./mean(Imag1_0p8uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p4uT(ii)=nansum(nansum(double(nImag1_0p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_0p8uT(ii)=nansum(nansum(double(nImag1_0p8uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p4uT_fatsuppression=(Zspectra_0p4uT(:)-Zspectra_0p4uT(26))./(1-Zspectra_0p4uT(26))
%  Zspectra_0p8uT_fatsuppression=(Zspectra_0p8uT(:)-Zspectra_0p8uT(26))./(1-Zspectra_0p8uT(26))
 Zspectra_0p4uT_fatsuppression_GM=Zspectra_0p4uT;
 Zspectra_0p8uT_fatsuppression_GM=Zspectra_0p8uT;

Zspectra_DPS_0p4uT_fatsuppression_GM=1./(1+(1./Zspectra_0p8uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p4uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p4uT_fatsuppression_GM-1./Zspectra_DPS_0p4uT_fatsuppression_GM;
MTRDSP_0p4uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p4uT_fatsuppression_GM-Zspectra_0p4uT_fatsuppression_GM;









%0p5uT and 1p0uT
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
    Zspectra_0p5uT(ii)=nansum(nansum(double(nImag1_0p5uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1uT(ii)=nansum(nansum(double(nImag1_1uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p5uT_fatsuppression=(Zspectra_0p5uT(:)-Zspectra_0p5uT(26))./(1-Zspectra_0p5uT(26))
%  Zspectra_1uT_fatsuppression=(Zspectra_1uT(:)-Zspectra_1uT(26))./(1-Zspectra_1uT(26))
 Zspectra_0p5uT_fatsuppression_GM=Zspectra_0p5uT;
 Zspectra_1uT_fatsuppression_GM=Zspectra_1uT;

Zspectra_DPS_0p5uT_fatsuppression_GM=1./(1+(1./Zspectra_1uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p5uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p5uT_fatsuppression_GM-1./Zspectra_DPS_0p5uT_fatsuppression_GM;
MTRDSP_0p5uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p5uT_fatsuppression_GM-Zspectra_0p5uT_fatsuppression_GM;






%0p6uT and 1p2uT
ls=dir('Sub6\dcm_format\MR_901_CEST_0_6uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_901_CEST_0_6uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p6uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p6uT(:,:,ii)=double(Imag1_0p6uT(:,:,ii))./mean(Imag1_0p6uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub6\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_1301_CEST_1_2uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p2uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p2uT(:,:,ii)=double(Imag1_1p2uT(:,:,ii))./mean(Imag1_1p2uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p6uT(ii)=nansum(nansum(double(nImag1_0p6uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p2uT(ii)=nansum(nansum(double(nImag1_1p2uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p6uT_fatsuppression_GM=Zspectra_0p6uT;
 Zspectra_1p2uT_fatsuppression_GM=Zspectra_1p2uT;

Zspectra_DPS_0p6uT_fatsuppression_GM=1./(1+(1./Zspectra_1p2uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p6uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p6uT_fatsuppression_GM-1./Zspectra_DPS_0p6uT_fatsuppression_GM;
MTRDSP_0p6uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p6uT_fatsuppression_GM-Zspectra_0p6uT_fatsuppression_GM;





%0p7 and 1p4
ls=dir('Sub6\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_1001_CEST_0_7uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_0p7uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_0p7uT(:,:,ii)=double(Imag1_0p7uT(:,:,ii))./mean(Imag1_0p7uT(:,:,[1,2,50,51]),3);
end

ls=dir('Sub6\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points');
    mniDir = ['Sub6\dcm_format\MR_1401_CEST_1_4uT_3_6s_51points\'];
for ij = 3:length(ls)
Imag1_1p4uT(:,:,ij-2)=dicomread([mniDir, ls(ij).name]);
end
for ii=1:51
nImag1_1p4uT(:,:,ii)=double(Imag1_1p4uT(:,:,ii))./mean(Imag1_1p4uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
    Zspectra_0p7uT(ii)=nansum(nansum(double(nImag1_0p7uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
    Zspectra_1p4uT(ii)=nansum(nansum(double(nImag1_1p4uT(:,:,ii)).*roi_brain_GM))./nansum(nansum(roi_brain_GM));
end

%  Zspectra_0p6uT_fatsuppression=(Zspectra_0p6uT(:)-Zspectra_0p6uT(26))./(1-Zspectra_0p6uT(26))
%  Zspectra_1p2uT_fatsuppression=(Zspectra_1p2uT(:)-Zspectra_1p2uT(26))./(1-Zspectra_1p2uT(26))
 Zspectra_0p7uT_fatsuppression_GM=Zspectra_0p7uT;
 Zspectra_1p4uT_fatsuppression_GM=Zspectra_1p4uT;

Zspectra_DPS_0p7uT_fatsuppression_GM=1./(1+(1./Zspectra_1p4uT_fatsuppression_GM-1).*((wl)/(wh))^2);
AREXDSP_0p7uT_fatsuppression_GM(:,rat_num)=1./Zspectra_0p7uT_fatsuppression_GM-1./Zspectra_DPS_0p7uT_fatsuppression_GM;
MTRDSP_0p7uT_fatsuppression_GM(:,rat_num)=Zspectra_DPS_0p7uT_fatsuppression_GM-Zspectra_0p7uT_fatsuppression_GM;


%Fig. S12a
% optimizing
sub_range=[1:6];
figure (1)
hold on
plot(flip(nanmean((MTRDSP_0p3uT_fatsuppression_WM(:,sub_range)-MTRDSP_0p3uT_fatsuppression_WM(46,sub_range)),2)));
plot(flip(nanmean((MTRDSP_0p4uT_fatsuppression_WM(:,sub_range)-MTRDSP_0p4uT_fatsuppression_WM(46,sub_range)),2)));
plot(flip(nanmean((MTRDSP_0p5uT_fatsuppression_WM(:,sub_range)-MTRDSP_0p5uT_fatsuppression_WM(46,sub_range)),2)));
plot(flip(nanmean((MTRDSP_0p6uT_fatsuppression_WM(:,sub_range)-MTRDSP_0p6uT_fatsuppression_WM(46,sub_range)),2)));
plot(flip(nanmean((MTRDSP_0p7uT_fatsuppression_WM(:,sub_range)-MTRDSP_0p7uT_fatsuppression_WM(46,sub_range)),2)));


%Fig. S12a
% optimizing
sub_range=[1:6];
figure (2)
hold on
plot(flip(nanmean((MTRDSP_0p3uT_fatsuppression_GM(:,sub_range)-MTRDSP_0p3uT_fatsuppression_GM(46,sub_range)),2)));
plot(flip(nanmean((MTRDSP_0p4uT_fatsuppression_GM(:,sub_range)-MTRDSP_0p4uT_fatsuppression_GM(46,sub_range)),2)));
plot(flip(nanmean((MTRDSP_0p5uT_fatsuppression_GM(:,sub_range)-MTRDSP_0p5uT_fatsuppression_GM(46,sub_range)),2)));
plot(flip(nanmean((MTRDSP_0p6uT_fatsuppression_GM(:,sub_range)-MTRDSP_0p6uT_fatsuppression_GM(46,sub_range)),2)));
plot(flip(nanmean((MTRDSP_0p7uT_fatsuppression_GM(:,sub_range)-MTRDSP_0p7uT_fatsuppression_GM(46,sub_range)),2)));


