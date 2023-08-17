function f=lambda(C1,C2,C3,C4,T,Tc) 
    f =10^-3*C1.*(1-T/Tc).^(C2+C3.*T/Tc+C4.*(T/Tc).^2) ;
end
