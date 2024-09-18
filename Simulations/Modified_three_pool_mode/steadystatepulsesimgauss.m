function p = steadystatepulsesimgauss(dw, ksw1,ksw2, ksw3,ksw4, ksw5, kmw, mnots1,mnots2, mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S1,R2S2,R2S3,R2S4,R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4,sep5, duration, curve, angle, TR, time, exciteflag, excitewait, exciteduration, exciteangle, spoil)
p = zeros(time,19);
init = [0;0;mnots1;0;0;mnotw;mnotm;0;0;mnots2; 0; 0; mnots3; 0; 0; mnots4; 0; 0; mnots5];
t = 1;
steadycheck = 100000;
steadystate = 1;
while (steadystate == 1 && t < time)
   blah = pulsesim1(dw, ksw1,ksw2, ksw3,ksw4,ksw5, kmw, mnots1,mnots2,mnots3, mnots4, mnots5, mnotw, mnotm, R1S, R2S1,R2S2,R2S3,R2S4,R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4, sep5, duration, curve, angle, init, TR, exciteflag, excitewait, exciteduration, exciteangle); 
   %blah;
   p(t,:) = blah;
   if(spoil == 1)

   end
   init = blah;

   t = t+1;

end
p(t:time, :) = [];
t
end