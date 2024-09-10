  clear;
  
% pulse shape
gauss =[
    0
    0
    0
     13.999  
     14.990  
     16.042  
     17.158  
     18.343  
     19.598  
     20.927  
     22.335  
     23.824  
     25.399  
     27.062  
     28.819  
     30.673  
     32.628  
     34.689  
     36.860  
     39.145  
     41.549  
     44.077  
     46.732  
     49.521  
     52.446  
     55.515  
     58.730  
     62.098  
     65.622  
     69.309  
     73.162  
     77.187  
     81.389  
     85.773  
     90.343  
     95.104  
    100.062  
    105.220  
    110.583  
    116.155  
    121.942  
    127.946  
    134.172  
    140.624  
    147.305  
    154.219  
    161.369  
    168.757  
    176.387  
    184.261  
    192.380  
    200.747  
    209.363  
    218.228  
    227.344  
    236.711  
    246.329  
    256.196  
    266.312  
    276.676  
    287.285  
    298.137  
    309.229  
    320.558  
    332.119  
    343.909  
    355.921  
    368.150  
    380.591  
    393.236  
    406.078  
    419.109  
    432.320  
    445.703  
    459.248  
    472.945  
    486.783  
    500.751  
    514.836  
    529.028  
    543.311  
    557.674  
    572.102  
    586.582  
    601.097  
    615.633  
    630.175  
    644.706  
    659.210  
    673.670  
    688.068  
    702.389  
    716.614  
    730.726  
    744.706  
    758.537  
    772.201  
    785.679  
    798.953  
    812.006  
    824.818  
    837.372  
    849.651  
    861.636  
    873.311  
    884.658  
    895.659  
    906.300  
    916.564  
    926.434  
    935.897  
    944.937  
    953.541  
    961.694  
    969.385  
    976.600  
    983.329  
    989.561  
    995.285  
   1000.492  
   1005.175  
   1009.324  
   1012.935  
   1016.000  
   1018.514  
   1020.474  
   1021.877  
   1022.719  
   1023.000  
   1022.719  
   1021.877  
   1020.474  
   1018.514  
   1016.000  
   1012.935  
   1009.324  
   1005.175  
   1000.492  
    995.285  
    989.561  
    983.329  
    976.600  
    969.385  
    961.694  
    953.541  
    944.937  
    935.897  
    926.434  
    916.564  
    906.300  
    895.659  
    884.658  
    873.311  
    861.636  
    849.651  
    837.372  
    824.818  
    812.006  
    798.953  
    785.679  
    772.201  
    758.537  
    744.706  
    730.726  
    716.614  
    702.389  
    688.068  
    673.670  
    659.210  
    644.706  
    630.175  
    615.633  
    601.097  
    586.582  
    572.102  
    557.674  
    543.311  
    529.028  
    514.836  
    500.751  
    486.783  
    472.945  
    459.248  
    445.703  
    432.320  
    419.109  
    406.078  
    393.236  
    380.591  
    368.150  
    355.921  
    343.909  
    332.119  
    320.558  
    309.229  
    298.137  
    287.285  
    276.676  
    266.312  
    256.196  
    246.329  
    236.711  
    227.344  
    218.228  
    209.363  
    200.747  
    192.380  
    184.261  
    176.387  
    168.757  
    161.369  
    154.219  
    147.305  
    140.624  
    134.172  
    127.946  
    121.942  
    116.155  
    110.583  
    105.220  
    100.062  
     95.104  
     90.343  
     85.773  
     81.389  
     77.187  
     73.162  
     69.309  
     65.622  
     62.098  
     58.730  
     55.515  
     52.446  
     49.521  
     46.732  
     44.077  
     41.549  
     39.145  
     36.860  
     34.689  
     32.628  
     30.673  
     28.819  
     27.062  
     25.399  
     23.824  
     22.335  
     20.927  
     19.598  
     18.343  
     17.158  
     16.042  
     14.990  
     0
     0
     0
     ];

% setting p1 and p2 values for the approximating model of repeated rectangular pulses 
p1=0.416;
p2=0.295;

% set the low saturation power
 w1_field=0.5;

% set the pulse duration 50ms
 pulseduration=0.05;
% set duty cycle
 dutycycle=0.999;
% calculate duty cycle for the approximating model of repeated rectangular pulses 
 DC_hard=dutycycle*p1^2/p2;
% pulse TR
 TR=pulseduration/dutycycle;

 % number of satuation pulses. The total saturation time is 3.6s.
 numb=fix(3.6/(TR));

% set frequency offset for each pool
% sep1 is amide, sep2 is amine, sep3 is guan, sep4 is NOE(-1.6), sep5 is
% NOE(-3.5)
sep1=3.5*127.7;
sep2=3*127.7;
sep3=2*127.7
sep4=-1.6*127.7;
sep5=-3.5*127.7;

%set pool concentration
% fs1 is amide, fs2 is amine, fs3 is guan, fs4 is NOE(-1.6), sep5 is
% NOE(-3.5), fm is MT.
fs1=0.0012;
fs2=0.003;
fs3=0.0005;
fs4=0.003;
fs5=0.008;
fm=0.15;

% set exchange rate for each pool (s^{-1})
% ksw1 is amide, fs2 is amine, fs3 is guan, fs4 is NOE(-1.6), fs5 is
% NOE(-3.5), kmw is MT-water coupling rate.
ksw1=100;
ksw2=5000;
ksw3=300;
ksw4=50;
ksw5=20;
kmw=25;

% set relaxations for each pool
R1S=1/1;
R2S1=1/0.006;
R2S2=1/0.02;
R2S3=1/0.02;
R2S4=1/0.003;
R2S5=1/0.0015;

R1W=1/1.;
R2W=1/0.09;
R1M=1/1.;
R2M=1/0.00004;

exciteduration=0;
exciteangle=0;

% set Z-spectra frequency offset
k= [-10,-8, -6, -5:0.25:5, 6, 8, 10]*127.7;
k=k';

% Bloch simulations of low power and high power CEST signals
w1_L=w1_field/dutycycle;
satangle=w1_L*42.6*360*pulseduration;
a25mspulse = runsteadysimgauss(ksw1, ksw2,ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3,fs4, fs5, 1, fm, R1S, R2S1,R2S2,R2S3,R2S4,R2S5, R1W, R2W, R1M, R2M,sep1*2*pi,sep2*2*pi,sep3*2*pi,sep4*2*pi, sep5*2*pi, pulseduration, gauss, satangle, TR, numb, 1, .002, exciteduration, exciteangle, k*2*pi, 1);
SL=a25mspulse(:,6);

w1_H=2*w1_field/dutycycle;
satangle=w1_H*42.6*360*pulseduration;
a25mspulse = runsteadysimgauss(ksw1, ksw2,ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3,fs4, fs5, 1, fm, R1S, R2S1,R2S2,R2S3,R2S4,R2S5, R1W, R2W, R1M, R2M,sep1*2*pi,sep2*2*pi,sep3*2*pi,sep4*2*pi, sep5*2*pi, pulseduration, gauss, satangle, TR, numb, 1, .002, exciteduration, exciteangle, k*2*pi, 1);
SH=a25mspulse(:,6);

% calculate auxilary signal
SA=1./(1+(1./SH-1)*(w1_L./w1_H)^2)





figure (1)
hold on
plot(SL)
plot(SH)


figure (2)
hold on
plot(SL)
plot(SA)


% plot MTRDSP and the residual MT corrected MTRDSP
figure (3)
hold on
%MTRDSP
plot(SA-SL)
%residual MT corrected MTRDSP 
plot(SA-SL-(SA(4)-SL(4)))



% plot AREXDSP and the residual MT corrected AREXDSP
figure (4)
hold on
%AREXDSP
plot((1./SL-1./SA).*(R1W+fm*R1M))
%residual MT corrected AREXDSP 
plot((1./SL-1./SA-(1./SL(4)-1./SA(4))).*(R1W+fm*R1M))


