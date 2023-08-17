clear all
clc

%Datos posibles de variar -------------------------------------------------
Gs0 = 0.956341779735860 ;
N = [0.8 0.9 1.1 1.25]*100 ;
Gs_v = Gs0*N/100 ;

%Fracción de las pérdidas en cada zona
FQp1 = 0.15;
FQp2 = 0.65;
FQp3 = 0.20;

%Datos de las corrientes --------------------------------------------------
TW1 = 20; %°C
w1 = 20; %En masa %
X1 = w1/(100-w1); %Fracción base seca

TW2 = 120; %°C
w2 = 0.3; %En masa %
X2 = w2/(100-w2); %Fracción base seca
W2 = 450/3600; %kg/h - flujo total
Ws = W2/(1 + X2); %kg solido seco/h 

TG2 = 155; %°C
Y2 = 0.01; %Fracción base seca
V_max = 1.6 ; % m/s

for o = 1:4
    
Gs = Gs_v(o) ; %°C - Valor asignado
V = 1.341265538111851 ; % m/s
%Datos de componentes -----------------------------------------------------

%Cps
global R CpWs CpH2OL

R = 8.314; %kJ/kmolK
CpWs = 837/1000; %kJ/kgK - sólido seco
CpH2OL = 4.187 ; %kJ/kgK

Tref = 0; %°C - Temp de referencia

%Balances globales --------------------------------------------------------

Qpf = @(Gs) 0.12*HG(TG2, Y2)*Gs ;

BM = @(Y1) Gs*(Y2-Y1) - Ws*(X2-X1);
BE = @(TG1,Y1) Gs*(HG(TG2,Y2)-HG(TG1,Y1)) - Qpf(Gs)  - Ws*(hW(TW2,X2) - hW(TW1,X1));

TG1 = 50; Y1 = 0.03;

for i=1:10
    X = fsolve(@(X) [BM(X(2)) ; BE(X(1),X(2))],[TG1;Y1]) ;
    TG1 = X(1) ;
    Y1 = X(2) ;
end

Qp = Qpf(Gs); %kJ/kgAS


tao = 10;
while tao > 0.0001
    
    FQp10 = FQp1;
    FQp20 = FQp2;
    FQp30 = FQp3;
    
    %Cálculos zona 3 ------------------------------------------------------

    %No hay vaporización en la zona 3
    YD = Y2; XB = X2; 
    %Balances
    BE3 = @(TWB,TGD) Gs*(HG(TG2,Y2) - HG(TGD,YD)) - Ws*(hW(TW2,X2) - hW(TWB,XB)) - Qp*FQp3;
    TbulbH = @(TWB,TGD) TWB - BulboHumedo(TGD,YD);
    x = fsolve(@(x) [BE3(x(1),x(2)) ; TbulbH(x(1),x(2))], [60,TW2]);
    TWB = x(1);
    TGD = x(2);

    DeltaTG3 = Ws*(CpWs + XB*CpH2OL)*(TW2-TWB)/(Gs*(CpAir_av((TG2+TGD)/2) + YD*CpH2OV_av((TG2+TGD)/2)));

    DeltaTml3 = ((TG2-TW2) - (TGD-TWB))/log((TG2-TW2)/(TGD-TWB));

    Ntog3 = DeltaTG3/DeltaTml3;


    %Cálculos zona 1 ------------------------------------------------------
    
    %No hay evaporación en la zona 1
    TWA = TWB; XA = X1; YC = Y1;
    BE1 = @(TGC) Gs*(HG(TGC,YC) - HG(TG1,Y1)) - Ws*(hW(TWA,XA) - hW(TW1,X1)) - Qp*FQp1;
    TGC = fsolve(BE1, TW1);

    DeltaTG1 = Ws*(CpWs + XA*CpH2OL)*(TWA-TW1)/(Gs*(CpAir_av((TG1+TGC)/2) + YC*CpH2OV_av((TG1+TGC)/2)));

    DeltaTml1 = ((TGC-TWA) - (TG1-TW1))/log((TGC-TWA)/(TG1-TW1));

    Ntog1 = DeltaTG1/DeltaTml1;


    %Cálculos zona 2 ------------------------------------------------------

    %Balances
    BE2 = @(TGC2) Gs*(CpAir_av((TGC2+TGD)/2) + CpH2OV_av((TGC2+TGD)/2)*(YD+YC)/2)*(TGD-TGC2)...
        - Ws*lambda(TWB)*(XA-XB) - Qp*FQp2;
    TGC2 = fsolve(@(TGC2) BE2(TGC2), TGC);


    DeltaTG2 = Ws*lambda(TWB)*(XA-XB)/(Gs*(CpAir_av((TGC+TGD)/2) + CpH2OV_av((TGC+TGD)/2)*(YD+YC)/2));

    DeltaTml2 = ((TGD-TWB) - (TGC-TWA))/log((TGD-TWB)/(TGC-TWA));

    Ntog2 = DeltaTG2/DeltaTml2;

    Ntog = Ntog1 + Ntog2 + Ntog3;
    FQp1 = Ntog1/Ntog;
    FQp2 = Ntog2/Ntog;
    FQp3 = Ntog3/Ntog;
    
    tao = (((FQp1 - FQp10)/FQp1)^2) + (((FQp2 - FQp20)/FQp2)^2) + (((FQp3 - FQp30)/FQp3)^2);
end


%Longitud unidades de transferencia ---------------------------------------

P1=101.325; %kPa
Yprom = (Y1 + Y2)/2;
vH = (R*(TG1+TG2+273.15*2)/(P1*2))*((1/29)+((Yprom)/(18))); % m3/kgAS

%Vmax = 1.6 m/s
G_A = V/vH ;
Gs_A = G_A/(1+Yprom) ;
A = Gs/Gs_A ;
Td = sqrt(A/pi)*2 ;


Gs_A = Gs/A; %Gasto kgAS/m2h



Ua = 237*(G_A^0.67)/(1000*Td); %kW/m2h

CsIn = (CpAir(TG2) + CpH2OV(TG2)*Y2);
CsOut = (CpAir(TG1) + CpH2OV(TG1)*Y1);
CsProm = (CsIn + CsOut)/2; %kJ/kgAS*K

Htog = Gs_A*CsProm/Ua; %m

Z = Ntog*Htog ;

Z_v(o)= Z;
TWA_v(o) = TWA;
TGC_v(o) = TGC;
TGD_v(o) = TGD;
Y1_v(o) = Y1;
TG1_v(o) = TG1;
Td_v(o) = Td ;
z(o,:) = [0,Ntog1*Htog, (Ntog2+Ntog1)*Htog, Ntog*Htog] ;
end
clc


%% Gráficas T vs Z
TW1_v = TW1*ones(size(TG1_v)) ;
TW2_v = TW2*ones(size(TG1_v)) ;
TG2_v = TG2*ones(size(TG1_v)) ;
TG = [TG1_v; TGC_v; TGD_v; TG2_v]' ;
TS = [TW1_v ; TWA_v ; TWA_v ; TW2_v]' ;

figure('Color','White')
t = tiledlayout(2,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

for i = 1:4
nexttile
hold on
plot(z(i,:),TS(i,:),'b-o',z(i,:),TG(i,:),'r-o')
title(sprintf("%.3f Gs",N(i)/100))
legend('Sólido','Gas','Location','NorthWest')
xlabel('Z, m')
ylabel('T, °C')
grid minor
end

%% Gráficas 

figure('Color','White')
t = tiledlayout(1,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(N,TG1_v,'-ok')
ylabel('TGS1, °C')
xlabel('%Gs')
grid minor

nexttile
plot(N,Y1_v,'-ok')
ylabel('Y1, °C')
xlabel('%Gs')
grid minor

figure('Color','White')
t = tiledlayout(1,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(N,Z_v,'-ok')
ylabel('Z, m')
xlabel('%Gs')
grid minor

nexttile
plot(N,Td_v,'-ok')
ylabel('D_T, m')
xlabel('%Gs')
grid minor

figure('Color','White')
t = tiledlayout(1,3) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(N,TWA_v,'-ok')
xlabel('%Gs')
ylabel('T_{WA} = T_{WB}, °C')
grid minor

nexttile
plot(N,TGC_v,'-ok')
xlabel('%Gs')
ylabel('T_{GC}, °C')
grid minor

nexttile
plot(N,TGD_v,'-ok')
xlabel('%Gs')
ylabel('T_{GD}, °C')
grid minor


%Funciones ----------------------------------------------------------------

function H = HG(T, Y)
    Tref = 0;
    H = (CpAir_av(T) + CpH2OV_av(T)*Y)*(T-Tref)+ lambda(Tref)*Y; %kJ/kgAS
end

function h = hW(T, X)
    global CpWs CpH2OL;
    Tref = 0;
    h = (CpWs + X*CpH2OL)*(T-Tref); %kJ/kgSS
end

%Temperatura de bulbo húmedo
function Tw = BulboHumedo(TG, Y)
    Cs = (CpAir(TG) + CpH2OV(TG)*Y);
    
    Obj1 = @(Tw) TG - (lambda(Tw)*(Ysat(Tw)-Y))/(Cs) - Tw;
    Tw = fsolve(Obj1, 30); %°C
end

function CpAguaProm = CpH2OV_av(T)
    Tref = 0; %°C
    CpAguaProm = integral(@(T) CpH2OV(T),Tref+273.15,T+273.15)/(T-Tref); %kJ/kgK
end

function CpAirProm = CpAir_av(T)
    Tref = 0; %°C
    CpAirProm = integral(@(T) CpAir(T),Tref+273.15,T+273.15)/(T-Tref); %kJ/kgK
end

function CpAire = CpAir(T)
    global R;
    CpAire = (3.355 + 0.575e-3*T - 0.016e5*T.^(-2))*R/29; %kJ/kgK
end

function CpAgua = CpH2OV(T)
    global R;
    CpAgua = (3.470 + 1.450e-3*T + 0.121e5*T.^(-2))*R/18; %kJ/kgK
end

%Humedad absoluta en la condición de bulbo húmedo
function ysat = Ysat(Tw)
    P = 760; %mmHg
    ysat = (18/29)*Pvap(Tw)/(P - Pvap(Tw));
end

%Calor de vaporización 
function f=lambda(T) 
    C1 = 5.2053*1e7;
    C2 = 0.3199;
    C3 = -0.212;
    C4 = 0.25795;
    Tr = (T+273.15)/647.14;
    f = C1*((1-Tr)^(C2+C3*Tr+C4*Tr*Tr))/18000; %kJ/kg
end

%Presión de vapor
function f = Pvap(T)  %Función Antoine
    A = 8.07131;
    B = 1730.630;
    C = 233.426;
    f = 10.^(A - B./(T + C)) ; %T en °C y P en mmHg
end
