function  p = pulsesim1(dw, ksw1,ksw2,ksw3, ksw4,ksw5, kmw, mnots1,mnots2, mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S1,R2S2,R2S3,R2S4,R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4,sep5, duration, curve, angle, init, TR, exciteflag, excitewait, exciteduration, exciteangle)




w1 = 0;
init(1)=0;
init(2)=0;
init(4)=0;
init(5)=0;
init(8)=0;
init(9)=0;
init(11)=0;
init(12)=0;
init(14)=0;
init(15)=0;
init(17)=0;
init(18)=0;


 [~,blah] = pulsesolv1(w1, dw, ksw1, ksw2,ksw3,ksw4,ksw5, kmw, mnots1,mnots2,mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S1,R2S2,R2S3,R2S4,R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4, sep5, init, TR-duration);

index = size(blah);
init = blah(index(1), :);   


init(1)=0;
init(2)=0;
init(4)=0;
init(5)=0;
init(8)=0;
init(9)=0;
init(11)=0;
init(12)=0;
init(14)=0;
init(15)=0;
init(17)=0;
init(18)=0;


for t=1:1:256

    w1=getsatpulse(curve, angle, t, duration, TR);
   [p,blah] = pulsesolv1(w1, dw, ksw1, ksw2,ksw3,ksw4,ksw5, kmw, mnots1,mnots2,mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S1,R2S2,R2S3,R2S4,R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4, sep5, init, duration/256);

   index = size(blah); 
   init = blah(index(1), :);   
end



p = blah(index(1), :);
end