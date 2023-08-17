function exceso = HE(Tx, x1) %Función para determinar la entalpía en exceso
    T_K = Tx + 273.15 ; %K
    a0 =-3.63868e5 + 1.83829e3.*T_K - 2.32763.*T_K.^2 ; 
    a1 = 9.25982e5 - 4.83586e3.*T_K + 6.37228.*T_K.^2 ;
    a2 = -14.04894e5 + 7.51661e3.*T_K - 10.11280.*T_K.^2 ;
    a3 = 10.91318e5 - 5.89498e3.*T_K + 7.98868.*T_K.^2 ;
    a4 = -2.79986e5 + 1.50557e3.*T_K - 2.03127.*T_K.^2  ;
    exceso = x1.*(1-x1).*(a0 + a1.*x1.^0.5+a2.*x1.^1.5+a3.*x1.^2.5+a4.*x1.^4.5) ; %J/mol
end