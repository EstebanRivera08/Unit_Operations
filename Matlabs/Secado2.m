clear
clc

%Información del sólido
w1 = 0.08 ; w2 = 0.005 ; %kgH2O/kgMix
TS1 = 27 ; TS2 = 150 ; %°C
dp = 200 ; %mu m
rho = 1300 ; %kg/m3
CpSS = 0.837 ; %J/kg K
mS2 = 0.63 ; %kg/s

%Información del Gas
G_Amax = 0.70 ; %kg/m2
yCO2 = 2.5/100 ; yO2 = 14.7/100 ; yN2 = 76/100 ; yH2O = 6.8/100 ;
y2 = [yCO2, yO2, yN2] ;
TG2 = 480 ; %°C

%Propiedades
%Gas
yGs = 1-yH2O ; %Fracción AS
R = 8.3145 ; %kJ/kgK
Tref = 0 ; %°C
HevapH2O = 2502.3 ; %kJ/kg

CpCO2 = @(T) 5.457 + 1.045e-3*T - 1.157e5*T.^(-2) ; 
CpO2 = @(T) 3.639 + 0.506e-3*T - 0.227e5*T.^(-2); 
CpN2 = @(T) 3.280 + 0.593e-3*T + 0.040e5*T.^(-2) ;
CpAir = @(T) 3.355 + 0.575e-3*T - 0.016e5*T.^(-2) ;
CpH2O = @(T) 3.470 + 1.450e-3*T + 0.121e5*T.^(-2) ;
Cp = @(T) [integral(CpCO2,Tref+273.15,T+273.15)/(T-Tref),...
    integral(CpO2,Tref+273.15,T+273.15)/(T-Tref),...
    integral(CpN2,Tref+273.15,T+273.15)/(T-Tref)]*R ; %kJ/kmolK
CpH2O = @(T) integral(CpH2O,Tref+273.15,T+273.15)/(T-Tref)*R ; %kJ/kmolK
CpH2OL = 4.187 ; %kJ/kgK

MWCO2 = 44 ; MWO2 = 32 ; MWN2 = 28 ; MWH2O = 18 ; %kg/kmol
MWGs = (yCO2*MWCO2 + yO2*MWO2 + yN2*MWN2)/yGs ; %kgAS/kmol1

%Balance de materia
X1 = w1/(1-w1) ; X2 = w2/(1-w2) ; %kgH2O/KgSs
Ss = mS2*(1-w2) ; %kgSs
Y2 = yH2O*MWH2O/(yGs*MWGs) ; %kgH2O/KgAS

BM = @(Gs,Y1) Gs*(Y2-Y1) - Ss*(X2-X1) ;

%Balance de energía
CpGs2 = (y2*Cp(TG2)')/(yGs*MWGs) ; CpGs1 = @(TG1) (y2*Cp(TG1)')/(yGs*MWGs) ;%kJ/kgAS K
CpH2O2 = CpH2O(TG2)/MWH2O ; CpH2O1 = @(TG1) CpH2O(TG1)/MWH2O ;%kJ/kgH2O K

HG2 = (CpGs2 + CpH2O2*Y2)*(TG2-Tref)+ HevapH2O*Y2 ; %kJ/kgAS
HG1 = @(TG1,Y1) (CpGs1(TG1) + CpH2O1(TG1)*Y1)*(TG1-Tref)+ HevapH2O*Y1 ; %kJ/kgAS
HS1 = (CpSS + X1*CpH2OL)*(TS1-Tref) ; %kJ/kgSS
HS2 = (CpSS + X2*CpH2OL)*(TS2-Tref) ; %kJ/kgSS
Q = @(Gs) 0.15*HG2*Gs ; %kJ/kgAS

BE = @(TG1,Gs,Y1) Gs*(HG2-HG1(TG1,Y1))-Q(Gs) - Ss*(HS2-HS1) ;

% Si TG1 = 120°C
TG1 = 120 ;

X = fsolve(@(X) [BM(X(1),X(2)) ; BE(120,X(1),X(2))],[100 ; 0.05]) ;
Gs = X(1) ;
Y1 = X(2) ;

ha = (Y1+Y2)/2 


