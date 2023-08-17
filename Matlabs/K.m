%Coeficiente de distribución - con Antoine
function f = K(i, T)
    global P;
    A = [6.80398 6.823 6.89385 7.2418];
    B = [803.81 1270.763 1264.37 1675.3];
    C = [246.99 221.416 216.636 200.01];
    
    f = (10^(A(i) - B(i)./(T + C(i))))*101.325/(760*P);
end