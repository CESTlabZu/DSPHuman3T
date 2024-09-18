clear all; close all; clc;
nx=144;
ny=144;
wl=0.5;
wh=1;

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


FreqArray = [-10,-8,-6,-5:0.25:5,6,8,10];
deta_cs =0.025;                                      
range_cs = -10:deta_cs:10;                             
interp1NOE = 801;                                  
referenceCriticalValue = 100;                     
cest1_0 = 0;         



% load ROI
load('data\roi_brain');
% read T1 map

sub_num=1;


for ii=1:51
nImag1_0p5uT(:,:,ii)=double(Imag1_0p5uT(:,:,ii))./mean(Imag1_0p5uT(:,:,[1,2,50,51]),3);
end

for ii=1:51
nImag1_1p0uT(:,:,ii)=double(Imag1_1p0uT(:,:,ii))./mean(Imag1_1p0uT(:,:,[1,2,50,51]),3);
end


% B0 correction 0.5uT
[Length,Width,Freq] = size(nImag1_0p5uT); 
nImag1_0p5uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_0p5uT(xx,yy,3:49).*roi_whole(xx,yy))';
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
[Length,Width,Freq] = size(nImag1_1p0uT); 
nImag1_1p0uT_corr(Length,Width,47)=0;
InterplCESTData(Length,Width,interp1NOE)=0;
for xx = 1:Length
    for yy = 1:Width
            Raw_Z = squeeze(nImag1_1p0uT(xx,yy,3:49).*roi_whole(xx,yy))';
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
                nImag1_1p0uT_corr(xx,yy,:) = TarValue_2;
            end
        end
end    

nImag1_DSP_0p5uT_corr=1./(1+(1./nImag1_1p0uT_corr-1).*((wl)/(wh))^2);
Image_MTR_APT_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,38)-nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_whole;
Image_MTR_NOE_0p5uT_corr(:,:,sub_num)=-(squeeze(nImag1_0p5uT_corr(:,:,10)-nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(nImag1_0p5uT_corr(:,:,44)-nImag1_DSP_0p5uT_corr(:,:,44))).*roi_whole;

Image_AREX_APT_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,38)-1./nImag1_DSP_0p5uT_corr(:,:,38))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_whole;
Image_AREX_NOE_0p5uT_corr(:,:,sub_num)=(squeeze(1./nImag1_0p5uT_corr(:,:,10)-1./nImag1_DSP_0p5uT_corr(:,:,10))-squeeze(1./nImag1_0p5uT_corr(:,:,44)-1./nImag1_DSP_0p5uT_corr(:,:,44))).*1000./double(Imag1_T1).*roi_whole;






for ii=1:47
Zspectra_0p5uT_eggwhite_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_eggwhite))./nansum(nansum(roi_eggwhite));
Zspectra_0p5uT_glu_pH7p2_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_glu_pH7p2))./nansum(nansum(roi_glu_pH7p2));
Zspectra_0p5uT_glu_pH7p0_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_glu_pH7p0))./nansum(nansum(roi_glu_pH7p0));
Zspectra_0p5uT_glu_pH6p5_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_glu_pH6p5))./nansum(nansum(roi_glu_pH6p5));
Zspectra_0p5uT_MnCl2_0p04mM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_MnCl2_0p04mM))./nansum(nansum(roi_MnCl2_0p04mM));
Zspectra_0p5uT_MnCl2_0p08mM_corr(ii,sub_num)=nansum(nansum(nImag1_0p5uT_corr(:,:,ii).*roi_MnCl2_0p08mM))./nansum(nansum(roi_MnCl2_0p08mM));

Zspectra_1p0uT_eggwhite_corr(ii,sub_num)=nansum(nansum(nImag1_1p0uT_corr(:,:,ii).*roi_eggwhite))./nansum(nansum(roi_eggwhite));
Zspectra_1p0uT_glu_pH7p2_corr(ii,sub_num)=nansum(nansum(nImag1_1p0uT_corr(:,:,ii).*roi_glu_pH7p2))./nansum(nansum(roi_glu_pH7p2));
Zspectra_1p0uT_glu_pH7p0_corr(ii,sub_num)=nansum(nansum(nImag1_1p0uT_corr(:,:,ii).*roi_glu_pH7p0))./nansum(nansum(roi_glu_pH7p0));
Zspectra_1p0uT_glu_pH6p5_corr(ii,sub_num)=nansum(nansum(nImag1_1p0uT_corr(:,:,ii).*roi_glu_pH6p5))./nansum(nansum(roi_glu_pH6p5));
Zspectra_1p0uT_MnCl2_0p04mM_corr(ii,sub_num)=nansum(nansum(nImag1_1p0uT_corr(:,:,ii).*roi_MnCl2_0p04mM))./nansum(nansum(roi_MnCl2_0p04mM));
Zspectra_1p0uT_MnCl2_0p08mM_corr(ii,sub_num)=nansum(nansum(nImag1_1p0uT_corr(:,:,ii).*roi_MnCl2_0p08mM))./nansum(nansum(roi_MnCl2_0p08mM));
end

 
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
  
 
threshold=4;
d_Zspectra_0p5uT_eggwhite_corr(4:44)=smoothdata(Zspectra_0p5uT_eggwhite_corr(4:44),"gaussian",threshold);
d_Zspectra_0p5uT_glu_pH7p2_corr(4:44)=smoothdata(Zspectra_0p5uT_glu_pH7p2_corr(4:44),"gaussian",threshold);
d_Zspectra_0p5uT_glu_pH7p0_corr(4:44)=smoothdata(Zspectra_0p5uT_glu_pH7p0_corr(4:44),"gaussian",threshold);
d_Zspectra_0p5uT_glu_pH6p5_corr(4:44)=smoothdata(Zspectra_0p5uT_glu_pH6p5_corr(4:44),"gaussian",threshold);
d_Zspectra_0p5uT_MnCl2_0p04mM_corr(4:44)=smoothdata(Zspectra_0p5uT_MnCl2_0p04mM_corr(4:44),"gaussian",threshold);
d_Zspectra_0p5uT_MnCl2_0p08mM_corr(4:44)=smoothdata(Zspectra_0p5uT_MnCl2_0p08mM_corr(4:44),"gaussian",threshold);

d_Zspectra_1p0uT_eggwhite_corr(4:44)=smoothdata(Zspectra_1p0uT_eggwhite_corr(4:44),"gaussian",threshold);
d_Zspectra_1p0uT_glu_pH7p2_corr(4:44)=smoothdata(Zspectra_1p0uT_glu_pH7p2_corr(4:44),"gaussian",threshold);
d_Zspectra_1p0uT_glu_pH7p0_corr(4:44)=smoothdata(Zspectra_1p0uT_glu_pH7p0_corr(4:44),"gaussian",threshold);
d_Zspectra_1p0uT_glu_pH6p5_corr(4:44)=smoothdata(Zspectra_1p0uT_glu_pH6p5_corr(4:44),"gaussian",threshold);
d_Zspectra_1p0uT_MnCl2_0p04mM_corr(4:44)=smoothdata(Zspectra_1p0uT_MnCl2_0p04mM_corr(4:44),"gaussian",threshold);
d_Zspectra_1p0uT_MnCl2_0p08mM_corr(4:44)=smoothdata(Zspectra_1p0uT_MnCl2_0p08mM_corr(4:44),"gaussian",threshold);

 
 
Zspectra_DSP_eggwhite_0p5uT_corr=1./(1+(1./d_Zspectra_1p0uT_eggwhite_corr-1).*((wl)/(wh))^2);
Zspectra_DSP_glu_pH7p2_0p5uT_corr=1./(1+(1./d_Zspectra_1p0uT_glu_pH7p2_corr-1).*((wl)/(wh))^2);
Zspectra_DSP_glu_pH7p0_0p5uT_corr=1./(1+(1./d_Zspectra_1p0uT_glu_pH7p0_corr-1).*((wl)/(wh))^2);
Zspectra_DSP_glu_pH6p5_0p5uT_corr=1./(1+(1./d_Zspectra_1p0uT_glu_pH6p5_corr-1).*((wl)/(wh))^2);
Zspectra_DSP_MnCl2_0p04mM_0p5uT_corr=1./(1+(1./d_Zspectra_1p0uT_MnCl2_0p04mM_corr-1).*((wl)/(wh))^2);
Zspectra_DSP_MnCl2_0p08mM_0p5uT_corr=1./(1+(1./d_Zspectra_1p0uT_MnCl2_0p08mM_corr-1).*((wl)/(wh))^2);


Zspectra_MTR_eggwhiteT_corr=Zspectra_DSP_eggwhite_0p5uT_corr-d_Zspectra_0p5uT_eggwhite_corr;
Zspectra_MTR_glu_pH7p2T_corr=Zspectra_DSP_glu_pH7p2_0p5uT_corr-d_Zspectra_0p5uT_glu_pH7p2_corr;
Zspectra_MTR_glu_pH7p0T_corr=Zspectra_DSP_glu_pH7p0_0p5uT_corr-d_Zspectra_0p5uT_glu_pH7p0_corr;
Zspectra_MTR_glu_pH6p5T_corr=Zspectra_DSP_glu_pH6p5_0p5uT_corr-d_Zspectra_0p5uT_glu_pH6p5_corr;
Zspectra_MTR_MnCl2_0p04mMT_corr=Zspectra_DSP_MnCl2_0p04mM_0p5uT_corr-d_Zspectra_0p5uT_MnCl2_0p04mM_corr;
Zspectra_MTR_MnCl2_0p08mMT_corr=Zspectra_DSP_MnCl2_0p08mM_0p5uT_corr-d_Zspectra_0p5uT_MnCl2_0p08mM_corr;

% R1
R1_map=1000./double(Imag1_T1);
R1_map(R1_map==inf)=nan;
R1_map(R1_map==-inf)=nan;
Value_R1_eggwhite=nansum(nansum(R1_map.*roi_eggwhite))./nansum(nansum(roi_eggwhite));
Value_R1_glu_pH7p2=nansum(nansum(R1_map.*roi_glu_pH7p2))./nansum(nansum(roi_glu_pH7p2));
Value_R1_glu_pH7p0=nansum(nansum(R1_map.*roi_glu_pH7p0))./nansum(nansum(roi_glu_pH7p0));
Value_R1_glu_pH6p5=nansum(nansum(R1_map.*roi_glu_pH6p5))./nansum(nansum(roi_glu_pH6p5));
Value_R1_MnCl2_0p04mM=nansum(nansum(R1_map.*roi_MnCl2_0p04mM))./nansum(nansum(roi_MnCl2_0p04mM));
Value_R1_MnCl2_0p08mM=nansum(nansum(R1_map.*roi_MnCl2_0p08mM))./nansum(nansum(roi_MnCl2_0p08mM));


Zspectra_AREX_eggwhiteT_corr=(1./d_Zspectra_0p5uT_eggwhite_corr-1./Zspectra_DSP_eggwhite_0p5uT_corr).*Value_R1_eggwhite;
Zspectra_AREX_glu_pH7p2T_corr=(1./d_Zspectra_0p5uT_glu_pH7p2_corr-1./Zspectra_DSP_glu_pH7p2_0p5uT_corr).*Value_R1_glu_pH7p2;
Zspectra_AREX_glu_pH7p0T_corr=(1./d_Zspectra_0p5uT_glu_pH7p0_corr-1./Zspectra_DSP_glu_pH7p0_0p5uT_corr).*Value_R1_glu_pH7p0;
Zspectra_AREX_glu_pH6p5T_corr=(1./d_Zspectra_0p5uT_glu_pH6p5_corr-1./Zspectra_DSP_glu_pH6p5_0p5uT_corr).*Value_R1_glu_pH6p5;
Zspectra_AREX_MnCl2_0p04mMT_corr=(1./d_Zspectra_0p5uT_MnCl2_0p04mM_corr-1./Zspectra_DSP_MnCl2_0p04mM_0p5uT_corr).*Value_R1_MnCl2_0p04mM;
Zspectra_AREX_MnCl2_0p08mMT_corr=(1./d_Zspectra_0p5uT_MnCl2_0p08mM_corr-1./Zspectra_DSP_MnCl2_0p08mM_0p5uT_corr).*Value_R1_MnCl2_0p08mM;

for ii=1:23
MTRasym_0p5uT_eggwhite_corr(ii)=d_Zspectra_0p5uT_eggwhite_corr(24-ii)-d_Zspectra_0p5uT_eggwhite_corr(24+ii);
MTRasym_0p5uT_glu_pH7p2_corr(ii)=d_Zspectra_0p5uT_glu_pH7p2_corr(24-ii)-d_Zspectra_0p5uT_glu_pH7p2_corr(24+ii);
MTRasym_0p5uT_glu_pH7p0_corr(ii)=d_Zspectra_0p5uT_glu_pH7p0_corr(24-ii)-d_Zspectra_0p5uT_glu_pH7p0_corr(24+ii);
MTRasym_0p5uT_glu_pH6p5_corr(ii)=d_Zspectra_0p5uT_glu_pH6p5_corr(24-ii)-d_Zspectra_0p5uT_glu_pH6p5_corr(24+ii);
MTRasym_0p5uT_MnCl2_0p04mM_corr(ii)=d_Zspectra_0p5uT_MnCl2_0p04mM_corr(24-ii)-d_Zspectra_0p5uT_MnCl2_0p04mM_corr(24+ii);
MTRasym_0p5uT_MnCl2_0p08mM_corr(ii)=d_Zspectra_0p5uT_MnCl2_0p08mM_corr(24-ii)-d_Zspectra_0p5uT_MnCl2_0p08mM_corr(24+ii);
end

fMTRasym_0p5uT_eggwhite_corr=flip(MTRasym_0p5uT_eggwhite_corr);
fMTRasym_0p5uT_glu_pH7p2_corr=flip(MTRasym_0p5uT_glu_pH7p2_corr);
fMTRasym_0p5uT_glu_pH7p0_corr=flip(MTRasym_0p5uT_glu_pH7p0_corr);
fMTRasym_0p5uT_glu_pH6p5_corr=flip(MTRasym_0p5uT_glu_pH6p5_corr);
fMTRasym_0p5uT_MnCl2_0p04mM_corr=flip(MTRasym_0p5uT_MnCl2_0p04mM_corr);
fMTRasym_0p5uT_MnCl2_0p08mM_corr=flip(MTRasym_0p5uT_MnCl2_0p08mM_corr);



figure (1)
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




figure (2)
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


figure (3)
hold on
plot(flip(Zspectra_MTR_eggwhiteT_corr))
plot(flip(Zspectra_MTR_glu_pH7p2T_corr))
plot(flip(Zspectra_MTR_glu_pH7p0T_corr))
plot(flip(Zspectra_MTR_glu_pH6p5T_corr))
plot(flip(Zspectra_MTR_MnCl2_0p04mMT_corr))
plot(flip(Zspectra_MTR_MnCl2_0p08mMT_corr))

plot(fMTRasym_0p5uT_glu_pH7p2_corr)
plot(fMTRasym_0p5uT_glu_pH7p0_corr)
plot(fMTRasym_0p5uT_glu_pH6p5_corr)
plot(flip(1-Zspectra_0p5uT_MnCl2_0p04mM_corr))
plot(flip(1-Zspectra_0p5uT_MnCl2_0p08mM_corr))



figure (4)
hold on
plot(flip(Zspectra_AREX_eggwhiteT_corr))
plot(flip(Zspectra_AREX_glu_pH7p2T_corr))
plot(flip(Zspectra_AREX_glu_pH7p0T_corr))
plot(flip(Zspectra_AREX_glu_pH6p5T_corr))
plot(flip(Zspectra_AREX_MnCl2_0p04mMT_corr))
plot(flip(Zspectra_AREX_MnCl2_0p08mMT_corr))
