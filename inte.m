function [disp_t,mom_t,H,Uinf,Delta]= inte(Y,WSmean)

Uinf=mean(WSmean(end-5:1:end));
WS_norm=WSmean/Uinf;

Ycor=Y;
for i=1:length(WS_norm)
    if WS_norm(i)>0.99
        j=i;
        break;
    end
end

WS_norm=transpose(WS_norm);
Ycor=transpose(Ycor);
Delta=pchip(WS_norm(1:1:j),Ycor(1:1:j),0.99);
disp_t=trapz([0,Ycor(1:1:j)],(1-[0,WS_norm(1:1:j)]));
mom_t=trapz([0,Ycor(1:1:j+2)],[0,WS_norm(1:1:j+2)].*(1-[0,WS_norm(1:1:j+2)]));
H=disp_t/mom_t;
