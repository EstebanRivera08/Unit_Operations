clear all
clc

%Datos posibles de variar -------------------------------------------------

%Fracci�n de las p�rdidas en cada zona
FQp1 = 0.15;
FQp2 = 0.65;
FQp3 = 0.20;

%Datos de componentes -----------------------------------------------------

%Cps
global R;
R = 8.314; %kJ/kmolK

global CpWs CpH2OL;
CpWs = 837/1000; %kJ/kgK - s�lido seco
CpH2OL = 4.187 ; %kJ/kgK

Tref = 0; %�C - Temp de referencia

n = 10 ;
Z_v = zeros(n,1);
Td_v = Z_v;
Gs_v = Z_v;
TWA_v = Z_v;
TGC_v = Z_v;
TGD_v = Z_v;
Y1_v = Z_v;
TG1_v = Z_v;


%Variaci�n de TG1
for o = 0:n-1

%Datos de las corrientes --------------------------------------------------
TW1 = 20; %�C
w1 = 20; %En masa %
X1 = w1/(100-w1); %Fracci�n base seca

TW2 = 120; %�C
w2 = 0.3; %En masa %
X2 = w2/(100-w2); %Fracci�n base seca
W2 = 450/3600; %kg/h - flujo total
Ws = W2/(1 + X2); %kg solido seco/h 

TG2 = 155; %�C
Y2 = 0.01; %Fracci�n base seca

TG1 = 40+o*8; %�C - Valor asignado


%Balances globales --------------------------------------------------------

Qpf = @(Gs) 0.12*HG(TG2, Y2)*Gs ;

BM = @(Gs,Y1) Gs*(Y2-Y1) - Ws*(X2-X1);
BE = @(Gs,Y1) Gs*(HG(TG2,Y2)-HG(TG1,Y1)) - Qpf(Gs)  - Ws*(hW(TW2,X2) - hW(TW1,X1));

Gs = 1; Y1 = 0.03;

for i=1:5
    X = fsolve(@(X) [BM(X(1),X(2)) ; BE(X(1),X(2))],[Gs;Y1]) ;
    Gs = X(1) ;
    Y1 = X(2) ;
end

Qp = Qpf(Gs);

tao = 10;
while tao > 0.0001
    
    FQp10 = FQp1;
    FQp20 = FQp2;
    FQp30 = FQp3;
    
    %C�lculos zona 3 ------------------------------------------------------

    %No hay vaporizaci�n en la zona 3
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


    %C�lculos zona 1 ------------------------------------------------------
    
    %No hay evaporaci�n en la zona 1
    TWA = TWB; XA = X1; YC = Y1;
    BE1 = @(TGC) Gs*(HG(TGC,YC) - HG(TG1,Y1)) - Ws*(hW(TWA,XA) - hW(TW1,X1)) - Qp*FQp1;
    TGC = fsolve(BE1, TW1);

    DeltaTG1 = Ws*(CpWs + XA*CpH2OL)*(TWA-TW1)/(Gs*(CpAir_av((TG1+TGC)/2) + YC*CpH2OV_av((TG1+TGC)/2)));

    DeltaTml1 = ((TGC-TWA) - (TG1-TW1))/log((TGC-TWA)/(TG1-TW1));

    Ntog1 = DeltaTG1/DeltaTml1;


    %C�lculos zona 2 ------------------------------------------------------

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
vH = (R*(TG1+TG2+273.15*2)/(P1*2))*((1/29)+(Yprom)/(18)); % m3/kgAS

%Vmax = 1.6 m/s
V = 1.3413; %m/s
A = Gs*vH/V; %m2
Td = sqrt(A/pi)*2;

Gs_A = Gs/A; %Gasto kgAS/m2h

G_A = Gs_A*(1 + Yprom); %kg/m2h

Ua = 237*(G_A^0.67)/(1000*Td); %kW/m2h

CsIn = (CpAir(TG2) + CpH2OV(TG2)*Y2);
CsOut = (CpAir(TG1) + CpH2OV(TG1)*Y1);
CsProm = (CsIn + CsOut)/2; %kJ/kgAS*K

Htog = Gs_A*CsProm/Ua; %m

Z = Ntog*Htog;

Z_v(o+1)= Z;
Td_v(o+1) = Td;
Gs_v(o+1) = Gs;
TWA_v(o+1) = TWA;
TGC_v(o+1) = TGC;
TGD_v(o+1) = TGD;
Y1_v(o+1) = Y1;
TG1_v(o+1) = TG1 ;

z(o+1,:) = [0,Ntog1*Htog, (Ntog2+Ntog1)*Htog, Ntog*Htog] ;

end
clc
%Gr�ficos -----------------------------------------------------------------
%%

figure('Color','White')
t = tiledlayout(1,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(TG1_v,Z_v,'k')
xlabel('T_G1, �C')
ylabel('Z, m')
grid minor

nexttile
plot(TG1_v,Td_v,'k')
xlabel('T_G1, �C')
ylabel('D_T, m')
grid minor

figure('Color','White')
t = tiledlayout(1,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(TG1_v,Gs_v,'k')
xlabel('T_{G1}, �C')
ylabel('Gs, kgAS/s')
grid minor

nexttile
plot(TG1_v,Y1_v,'k')
xlabel('T_{G1}, �C')
ylabel('Y_1, kgH_2O/kgAS')
grid minor

figure('Color','White')
t = tiledlayout(1,3) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(TG1_v,TWA_v,'k')
xlabel('T_{G1}, �C')
ylabel('T_{WA} = T_{WB}, �C')
grid minor

nexttile
plot(TG1_v,TGC_v,'k')
xlabel('T_{G1}, �C')
ylabel('T_{GC}, �C')
grid minor

nexttile
plot(TG1_v,TGD_v,'k')
xlabel('T_{G1}, �C')
ylabel('T_{GD}, �C')
grid minor
%%
figure('Color','White')
TW1_v = TW1*ones(size(TG1_v)) ;
TW2_v = TW2*ones(size(TG1_v)) ;
TG2_v = TG2*ones(size(TG1_v)) ;
TS = [TG1_v TGC_v TGD_v TG2_v] ;
TG = [TW1_v TWA_v TWA_v TW2_v] ;
for i = 1 :3: n
hold on
plot(z(i,:),TS(i,:),'-o')
plot(z(i,:),TG(i,:),'-o')
end
grid minor
%%
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

%Temperatura de bulbo h�medo
function Tw = BulboHumedo(TG, Y)
    Cs = (CpAir(TG) + CpH2OV(TG)*Y);
    
    Obj1 = @(Tw) TG - (lambda(Tw)*(Ysat(Tw)-Y))/(Cs) - Tw;
    Tw = fsolve(Obj1, 30); %�C
end

function CpAguaProm = CpH2OV_av(T)
    Tref = 0; %�C
    CpAguaProm = integral(@(T) CpH2OV(T),Tref+273.15,T+273.15)/(T-Tref); %kJ/kgK
end

function CpAirProm = CpAir_av(T)
    Tref = 0; %�C
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

%Humedad absoluta en la condici�n de bulbo h�medo
function ysat = Ysat(Tw)
    P = 760; %mmHg
    ysat = (18/29)*Pvap(Tw)/(P - Pvap(Tw));
end

%Calor de vaporizaci�n 
function f=lambda(T) 
    C1 = 5.2053*1e7;
    C2 = 0.3199;
    C3 = -0.212;
    C4 = 0.25795;
    Tr = (T+273.15)/647.14;
    f = C1*((1-Tr)^(C2+C3*Tr+C4*Tr*Tr))/18000; %kJ/kg
end

%Presi�n de vapor
function f = Pvap(T)  %Funci�n Antoine
    A = 8.07131;
    B = 1730.630;
    C = 233.426;
    f = 10.^(A - B./(T + C)) ; %T en �C y P en mmHg
end
