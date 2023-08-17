
function [f1,f2] = funcion(XS,NS) 

%1 - Propano
m1 = 44; %g/mol
%2 - Metilciclohexano
m2 = 98.186; %g/mol
%3 - Heptano
m3 = 100.2019; %g/mol
%4 - Anilina
m4 = 93.12; %g/mol

%OPERACIÓN 1 - DESTILACIÓN SIMPLE MULTICOMPONENTE -------------------------
%Vectores -----------------------------------------------------------------

zF1 = zeros(3, 1); 
y1 = zeros(3, 1); 
x1 = zeros(3, 1);

%Entradas -----------------------------------------------------------------

global P
F1 = 2000*453.59; %mol/h
P = (20 + 14.69)*6.89476; %kPa
%P = 101.325;
P=70;
TF1 = 30; %°C
zF1(1) = 0.2; zF1(2) = 0.35; zF1(3) = 0.45; %mol

%Flujos -------------------------------------------------------------------

%Relación de flujos (V/F) ---- W
Tenedor = @(J, T) (K(1,T) - 1).*zF1(1)./(J.*(K(1,T) - 1) + 1)...
    + (K(2,T) - 1).*zF1(2)./(J.*(K(2,T) - 1) + 1)...
    + (K(3,T) - 1).*zF1(3)./(J.*(K(3,T) - 1) + 1);

J = fsolve(@(J) Tenedor(J,TF1), 0.2);

%Composiciones de salida
for i = 1:3
    x1(i) = zF1(i)/(J.*(K(i,TF1) - 1) + 1);
    y1(i) = K(i,TF1)*zF1(i)/(J.*(K(i,TF1) - 1) + 1);
end


%OPERACIÓN 2 - EXTRACCIÓN LIQ-LIQ - REFLUJO DE EXTRACTO -------------------
%Trabajo con base libre de solvente 

%Equilibrio ---------------------------------------------------------------

%En el refinado
R_Metilciclohexano = [0 9.2 18.6 22 33.8 40.9 46 59 67.2 71.6 ...
    73.6 83.3 88.1]';
R_Heptano = [92.9 83.1 73.4 69.8 57.6 50.4 45 30.7 22.8 ...
    18.2 16 5.4 0]';
R_Anilina = 100-R_Metilciclohexano-R_Heptano ;
%En el extracto
E_Metilciclohexano = [0 0.8 2.7 3 4.6 6 7.4 9.2 11.3 12.7 13.1 ...
    15.6 16.9]';
E_Heptano = [6.2 6 5.3 5.1 4.5 4 3.6 2.8 2.1 1.6 1.4 ...
    0.6 0]';
E_Anilina = 100-E_Metilciclohexano-E_Heptano ;

%Base libre
N_E = E_Anilina./(100-E_Anilina) ;
N_R = R_Anilina./(100-R_Anilina) ;
X = R_Metilciclohexano./(100-R_Anilina) ;
Y = E_Metilciclohexano./(100-E_Anilina) ;

%Funciones de interpolación
Y_func = @(x) interp1(X, Y, x); %Función para Y en términos de X 
X_func = @(y) interp1(Y, X, y); %Función para X en términos de Y 
N_R_func = @(x) interp1(X, N_R, x); %Función para N_R en términos de X
N_E_func = @(y) interp1(Y, N_E, y); %Función para N_E en términos de Y

%Variables ----------------------------------------------------------------

%Alimento
F_2 = F1*(1-J); %mol/h
%Alimento - fracciones en masa (Sin propano)
mF_mass = x1(2)*m2/(x1(2)*m2+x1(3)*m3);
hF_mass = x1(3)*m3/(x1(2)*m2+x1(3)*m3);
%Alimento - base libre
XF = mF_mass/(mF_mass+hF_mass);
NF = 0;
F2 = (x1(2)*m2+x1(3)*m3)*F_2/1000; %kg/h F - base libre de solvente (anilina) 

%Suposiciones -------------------------------------------------------------

%Fondo columna de destilación 2 - fracciones molares
hB_mol = 0.003; mB_mol = hB_mol*13.1/1.4; aB_mol = 1 - hB_mol - mB_mol; 
%Fondo columna de destilación 2 - fracciones en masa (Sin propano)
aB_mass = aB_mol*m4/(mB_mol*m2+hB_mol*m3+aB_mol*m4);
mB_mass = mB_mol*m2/(mB_mol*m2+hB_mol*m3+aB_mol*m4);
hB_mass = hB_mol*m3/(mB_mol*m2+hB_mol*m3+aB_mol*m4);
%Fondo columna de destilación 2 - base libre
XB = mB_mass/(mB_mass+hB_mass);
NB = aB_mass/(mB_mass+hB_mass);

%Destilado columna de destilación 2 - fracciones molares
hE_mol = 0.037; mE_mol = hE_mol*13.1/1.4; aE_mol = 1 - hE_mol - mE_mol; 
%Destilado columna de destilación 2 - fracciones en masa (Sin propano)
aE_mass = aE_mol*m4/(mE_mol*m2+hE_mol*m3+aE_mol*m4);
mE_mass = mE_mol*m2/(mE_mol*m2+hE_mol*m3+aE_mol*m4);
hE_mass = hE_mol*m3/(mE_mol*m2+hE_mol*m3+aE_mol*m4);
%Destilado columna de destilación 2 - base libre
%XE = mE_mass/(mE_mass+hE_mass);
XE = mE_mass/(mE_mass+hE_mass);
NE = aE_mass/(mE_mass+hE_mass);
XR0 = XE; NR0 = NE;
XC = XE; NC = NE;

%E1 es la suma de B y E, se usa palanca para determinarlo
YE1 = XE;
NE1 = N_E_func(YE1);

%Solvente(entrada a la columna) - base libre
% XS = 0.2;
% %Solvente(entrada a la columna) - fracciones en masa (Sin propano)
% aS_mass = 0.95;
% Obj1 = @(x) XS - (x/(1-aS_mass));
% mS_mass=fsolve(Obj1, 0.027);
% hS_mass=1-aS_mass-mS_mass;
% %Solvente(entrada a la columna) - base libre
% NS = aS_mass/(mS_mass+hS_mass);

%Refinado final - fracciones en masa (Sin propano)
aRNp_mass = 0.08;
mRNp_mass = 0.17;
hRNp_mass = 0.75;
%Refinado final - base libre
XRNp = mRNp_mass/(mRNp_mass+hRNp_mass);
NRNp = aRNp_mass/(mRNp_mass+hRNp_mass);

%Relación de reflujo mínima -----------------------------------------------

%Equilibrio que al extender toca al alimento
Obj2 = @(X) ((N_R_func(X) - N_E_func(Y_func(X)))/(X - Y_func(X)))...
    - ((N_R_func(X) - NF)/(X - XF));
XRmin = fsolve(Obj2, 0.4);

%Palanca para DeltaEmin
XDeltaEmin = XB;
Obj3 = @(X) ((Y_func(XRmin) - XDeltaEmin)/(N_E_func(Y_func(XRmin)) - X)) - ((XRmin - Y_func(XRmin))/(N_R_func(XRmin) - N_E_func(Y_func(XRmin))));
NDeltaEmin = fsolve(Obj3, 25);

%Doble palanca para DeltaRmin
Obj4 = @(X) [((XF - X(1))/(NF - X(2))) - ((XDeltaEmin - XF)/(NDeltaEmin - NF)),...
    ((NS - X(2))/(XS - X(1))) - ((NS - NRNp)/(XS - XRNp))];
DeltaRmin = fsolve(Obj4, [0.1, -20]);
XDeltaRmin = DeltaRmin(1);
NDeltaRmin = DeltaRmin(2);

RRmin = (NDeltaEmin - NE1)/(NE1 - NR0);

%Operatorias --------------------------------------------------------------

RR = 1.2*RRmin; %Relación de reflujo
NDeltaE = RR*(NE1 - NR0) + NE1;
XDeltaE = XDeltaEmin;

%Doble palanca para DeltaR
Obj6 = @(X) [((XF - X(1))/(NF - X(2))) - ((XDeltaE - XF)/(NDeltaE - NF)),...
    ((NS - X(2))/(XS - X(1))) - ((NS - NRNp)/(XS - XRNp))];
DeltaR = fsolve(Obj6, [XDeltaRmin - 0.1, NDeltaRmin - 10]);
XDeltaR = DeltaR(1);
NDeltaR = DeltaR(2);

%También se halla el X en el que se van a encontrar las líneas de operación
Obj7 = @(X) ((X - XF)/(N_R_func(X) - NF)) - ((XDeltaE - XF)/(NDeltaE - NF));
Xcorte = fsolve(Obj7, XRmin);

%Se determinan las parejas de puntos para la operatorias
X_opE = linspace(Xcorte, XR0-0.001, 20);
X_opR = linspace(XRNp+0.001, Xcorte, 20);
for i = 1 : 20
    Obj8_E = @(Y) ((XDeltaE - Y)./(NDeltaE - N_E_func(Y))) - ...
        ((XDeltaE - X_opE(i))./(NDeltaE - N_R_func(X_opE(i))));
    Y_opE(i) = fsolve(Obj8_E, 0.5);
    
    Obj8_R = @(Y) ((XDeltaR - Y)./(NDeltaR - N_E_func(Y))) - ...
        ((XDeltaR - X_opR(i))./(NDeltaR - N_R_func(X_opR(i))));
    Y_opR(i) = fsolve(Obj8_R, 0.1);
end

%Funciones de las lineas de operación
Op_R = @(X) interp1(X_opR, Y_opR, X, 'linear','extrap');
Op_E = @(X) interp1(X_opE, Y_opE, X, 'linear','extrap');

%Estimación del número de etapas y sus composiciones ----------------------

EtX(1) = XR0;
EtY(1) = Op_E(XR0);
EtY(2) = EtY(1);
EtX(2) = X_func(EtY(2));
i = 2;
Et = 1;
while EtX(i) > XRNp
    if EtX(i) > Xcorte
        EtX(i+1) = EtX(i);
        EtY(i+1) = Op_E(EtX(i+1));
        EtY(i+2) = EtY(i+1);
        EtX(i+2) = X_func(EtY(i+2));
    elseif EtX(i) > XRNp+0.001
        EtX(i+1) = EtX(i);
        EtY(i+1) = Op_R(EtX(i+1));
        EtY(i+2) = EtY(i+1);
        EtX(i+2) = X_func(EtY(i+2));
    else
        EtX(i+1) = EtX(i);
        EtY(i+1) = Op_R(EtX(i+1));
    end
    i = i + 2;
    Et = Et + 1;
end    

Et = size(EtX,2)/2 ; %Número de etapas ideales
  
for i=1:1:Et
    X_et(i) = EtX(2*i-1); %Para concentración en R0.....R10
    NR_et(i) = N_R_func(X_et(i));
    Y_et(i) = EtY(2*i-1); %Para concentración en E1.....E11  
    NE_et(i) = N_E_func(Y_et(i));
end  

X_et(Et+1) = XRNp;
NR_et(Et+1) = NRNp;
Y_et(Et+1) = XS;
NE_et(Et+1) = NS;

%Cálculo de flujos --------------------------------------------------------

DeltaR = F2*(XF - XDeltaE)/(XDeltaR - XDeltaE); %kg/h - base libre
DeltaE = F2 - DeltaR;

R_et(1) = DeltaE*RR; %Magnitud de R0 - base libre
E_et(1) = R_et(1)*(NDeltaE - NR_et(1))/(NDeltaE - NE_et(1)); %Magnitud de E1

for i=2:(Et)
    %Enriquecimiento
    if X_et(i) > Xcorte
        E_et(i) = DeltaE/(1-((NDeltaE - NE_et(i))/(NDeltaE - NR_et(i))));
        R_et(i) =  E_et(i) - DeltaE;   
    %Despojamiento
    else
        E_et(i) = DeltaR/(((NE_et(i) - NDeltaR)/(NR_et(i) - NDeltaR)) - 1);
        R_et(i) = E_et(i) + DeltaR;
    end
end

S = DeltaR/(((NDeltaR - NS)/(NDeltaR - NRNp)) - 1);
RNp = S + DeltaR;
E_et(Et+1) = S;
R_et(Et+1) = RNp;
E1 = E_et(1) ;
R0 = R_et(1) ;

%Balance en la columna de destilación -------------------------------------

%Balances en base libre
B = E1*(NE_et(1) - NE)/(NB - NE);
E = E1 - B;

C = E - R0;

%Balances del mezclador y el divisor de refinado --------------------------

%Balances en base libre
Z = S - B;
A = C*NC + (RNp - Z)*NRNp;
XSS = (B*XB + Z*XRNp)/S;
NSS = (B*NB + A + Z*NRNp)/S;

f1 = XSS-XS ;
f2 = NSS-NS ;

for i=1:Et
    B_et(i) = -E_et(i) + R_et(i) - R_et(i+1) + E_et(i+1);
end

end

%FUNCIONES ----------------------------------------------------------------


