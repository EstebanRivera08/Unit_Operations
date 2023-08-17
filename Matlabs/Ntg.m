clc
clear

yA1 = 0.042;
YA1 = yA1./(1-yA1) ;
abs = 0.99 ;
YA2 = (1-abs)*YA1  ;
yA2 = YA2/(1+YA2) ;

m = 10.7 ;
Fl = 1.41e-2 ;
Fg = 6.84e-4 ;
Ls_Gs = 16.521 ;
Gs_Ls = 1/Ls_Gs ;


k = @(yA) (yA./(1-yA) - YA2)*Gs_Ls ;

xA = @(yA) k(yA)./(1+k(yA)) ;

xAi = @(yAi) 1/m*yAi ;

yA = linspace(yA1,yA2,500) ;

Obj = @(yAi,yA) ((1-yAi)./(1-yA)) - ((1-xA(yA))./(1-xAi(yAi))).^(Fl/Fg) ;

for i = 1:size(yA,2)
    yAi(i) = fsolve(@(yAi) Obj(yAi,yA(i)),yA(i)) ;
end
clc

yAi_ = @(x) interp1(yA,yAi,x,'Linear','extrap') ;

f = @(yA) 1./((1-yA).*log((1-yAi_(yA))./(1-yA))) ;

Ntg = integral(f,yA2,yA1) 





