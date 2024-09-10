function [a,b] = pulsesolv1(w1, dw, ksw1,ksw2,ksw3,ksw4, ksw5,kmw, mnots1,mnots2,mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S1,R2S2,R2S3,R2S4,R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3, sep4,sep5, init, duration)
tspan = [0 duration];
y0 = init;

delter=2.3*127.7*2*pi;


ls=1./(pi*R2M)./(1+((dw-delter)./R2M).^2);
W=(w1)^2*pi*ls;  


[a,b] = ode45(@f, tspan, y0);
    function dydt = f(t,y)   
    dydt = [ ((-sep1 - dw)*y(2) - (R2S1 + ksw1)*y(1) + mnots1/mnotw*ksw1*y(4))
             ((sep1 + dw)*y(1) + w1*y(3) - (R2S1 + ksw1) * y(2) + mnots1/mnotw * ksw1 * y(5))
             (-w1*y(2) - R1S*(y(3) - mnots1) - ksw1*y(3) + mnots1/mnotw*ksw1*y(6))
             (-dw*y(5)-R2W*y(4)+ksw1*y(1)- mnots1/mnotw*ksw1*y(4)+ksw2*y(8)- mnots2/mnotw*ksw2*y(4))
             (dw*y(4) + 0*y(6)-R2W*y(5) + ksw1*y(2) - mnots1/mnotw*ksw1 *y(5)+ ksw2*y(9) - mnots2/mnotw*ksw2 *y(5)+ ksw3*y(12) - mnots3/mnotw*ksw3 *y(5)+ ksw4*y(15) - mnots4/mnotw*ksw4 *y(5)+ ksw5*y(18) - mnots5/mnotw*ksw5 *y(5))
             (-0*y(5) - R1W*(y(6) - mnotw) + ksw1*y(3) - mnots1/mnotw*ksw1*y(6)+ ksw2*y(10) - mnots2/mnotw*ksw2*y(6)+ ksw3*y(13) - mnots3/mnotw*ksw3*y(6)+ ksw4*y(16) - mnots4/mnotw*ksw4*y(6)+ ksw5*y(19) - mnots5/mnotw*ksw5*y(6)+kmw*y(7)-mnotm/mnotw*kmw*y(6))
             (- R1M*(y(7) - mnotm) - kmw*y(7) + mnotm/mnotw*kmw*y(6)-W*y(7))
             ((-sep2 - dw)*y(9) - (R2S2 + ksw2)*y(8) + mnots2/mnotw*ksw2*y(4))
             ((sep2 + dw)*y(8) + w1*y(10) - (R2S2 + ksw2) * y(9) + mnots2/mnotw * ksw2 * y(5))
             (-w1*y(9) - R1S*(y(10) - mnots2) - ksw2*y(10) + mnots2/mnotw*ksw2*y(6))    
             ((-sep3 - dw)*y(12) - (R2S3 + ksw3)*y(11) + mnots3/mnotw*ksw3*y(4))
             ((sep3 + dw)*y(11) + w1*y(13) - (R2S3 + ksw3) * y(12) + mnots3/mnotw * ksw3 * y(5))
             (-w1*y(12) - R1S*(y(13) - mnots3) - ksw3*y(13) + mnots3/mnotw*ksw3*y(6))   
             ((-sep4 - dw)*y(15) - (R2S4 + ksw4)*y(14) + mnots4/mnotw*ksw4*y(4))
             ((sep4 + dw)*y(14) + w1*y(16) - (R2S4 + ksw4) * y(15) + mnots4/mnotw * ksw4 * y(5))
             (-w1*y(15) - R1S*(y(16) - mnots4) - ksw4*y(16) + mnots4/mnotw*ksw4*y(6))   
             ((-sep5 - dw)*y(18) - (R2S5 + ksw5)*y(17) + mnots5/mnotw*ksw5*y(4))
             ((sep5 + dw)*y(17) + w1*y(19) - (R2S5 + ksw5) * y(18) + mnots5/mnotw * ksw5 * y(5))
             (-w1*y(18) - R1S*(y(19) - mnots5) - ksw5*y(19) + mnots5/mnotw*ksw5*y(6))   
               ];
             
           
    end
end
