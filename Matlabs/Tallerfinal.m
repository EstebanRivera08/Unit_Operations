clear all
clc

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
mB_mass = 0.02; hB_mass = mB_mass*1.4/13.1; aB_mass = 1 - hB_mass - mB_mass; 
%Fondo columna de destilación 2 - fracciones en masa (Sin propano)
aB_mol = (aB_mass/m4)/((mB_mass/m2)+(hB_mass/m3)+(aB_mass/m4));
mB_mol = (mB_mass/m2)/((mB_mass/m2)+(hB_mass/m3)+(aB_mass/m4));
hB_mol = (hB_mass/m3)/((mB_mass/m2)+(hB_mass/m3)+(aB_mass/m4));
%Fondo columna de destilación 2 - base libre
XB = mB_mass/(mB_mass+hB_mass);
NB = aB_mass/(mB_mass+hB_mass);

%Destilado columna de destilación 2 - fracciones molares
mE_mass = 0.90; hE_mass = mE_mass*1.4/13.1; aE_mass = 1 - hE_mass - mE_mass; 
%Destilado columna de destilación 2 - fracciones en masa (Sin propano)
aE_mol = (aE_mass/m4)/((mE_mass/m2)+(hE_mass/m3)+(aE_mass/m4));
mE_mol = (mE_mass/m2)/((mE_mass/m2)+(hE_mass/m3)+(aE_mass/m4));
hE_mol = (hE_mass/m3)/((mE_mass/m2)+(hE_mass/m3)+(aE_mass/m4));
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
XS = 0.12;
%Solvente(entrada a la columna) - fracciones en masa (Sin propano)
aS_mass = 0.95;
Obj1 = @(x) XS - (x/(1-aS_mass));
mS_mass=fsolve(Obj1, 0.027);
hS_mass=1-aS_mass-mS_mass;
%Solvente(entrada a la columna) - base libre
NS = aS_mass/(mS_mass+hS_mass);

%Refinado final - fracciones en masa (Sin propano)
aRNp_mass = 0.08;
mRNp_mass = 0.10;
hRNp_mass = 0.82;
%Refinado final - base libre
XRNp = mRNp_mass/(mRNp_mass+hRNp_mass);
NRNp = aRNp_mass/(mRNp_mass+hRNp_mass) ;

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

%Balance en la columna de destilación -------------------------------------

%Balances en base libre
B = E_et(1)*(NE_et(1) - NE)/(NB - NE);
E = E_et(1) - B;

C = E - R_et(1);

%Balances del mezclador y el divisor de refinado --------------------------

%Balances en base libre
Z = S - B;
A = C*NC + (RNp - Z)*NRNp;
Ve = RNp - Z;


for i=1:Et
    B_et(i) = -E_et(i) + R_et(i) - R_et(i+1) + E_et(i+1);
end



%% OPERACIÓN 3 - COLUMNA DE DESTILACIÓN -------------------------------------

%Datos --------------------------------------------------------------------

T_ref = 25; %C
Cpa_liq = 200; %kJ/kmolK
Cpa_gas = 149.06; %kJ/kmolK
Hvapa = 46.72*1000; %kJ/kmolK
Cpm_gas = 142.7; %kJ/kmolK
Cpm_liq = 184.84; %kJ/kmolK
Hvapm = 35.3*1000; %kJ/kmolK
Cph_gas = 191.5; %kJ/kmolK
Cph_liq = 224.721; %kJ/kmolK
Hvaph = 36.3*1000; %kJ/kmolK

m_h = 13.1/1.4; %Relación M-H

Cpmix_liq = (Cpm_liq*m_h+Cph_liq)/(m_h+1); %kJ/kmolK
Cpmix_gas = (Cpm_gas*m_h+Cph_gas)/(m_h+1); %kJ/kmolK
Hvapmix = (Hvapm*m_h+Hvaph)/(m_h+1); %kJ/kmolK

%Equilibrio ---------------------------------------------------------------

%Temperatura-composiciones
Obj9 = @(x,T) x*K(2,T)+(1-x)*K(4,T)-1;
x = 0:0.01:1;
for i=1:size(x,2)
    Tb(i) = fsolve(@(T) Obj9(x(i),T), 130) ;
    y(i) = K(2,Tb(i))*x(i) ;
end
y_func = @(u) interp1(x, y, u, 'linear','extrap'); %Función
x_func = @(u) interp1(y, x, u, 'linear','extrap'); %Función para x en términos de y 
T_func = @(u) interp1(x, Tb, u, 'linear','extrap'); %Función de T para cada x en equilibrio

%Entalpías
HL = (x.*Cpmix_liq.*(Tb-T_ref) + (1-x).*Cpa_liq.*(Tb-T_ref)) ; %kJ/kmol
HV = (y.*(Hvapmix + Cpmix_gas.*(Tb-T_ref)) + (1-y).*(Hvapa + Cpa_gas.*(Tb-T_ref))); %kJ/kmol
Fh = @(u) interp1(x,HL,u, 'linear','extrap');
FH = @(u) interp1(y,HV,u, 'linear','extrap');

%Entrada y salidas---------------------------------------------------------

P = 101.325; %kPa
aE1_mass = NE_et(1)/(1 + NE_et(1));
mE1_mass = YE1*(1-aE1_mass);
hE1_mass = 1 - aE1_mass - mE1_mass;
F3_mass = E_et(1)*(1 + NE_et(1)); %kg/h Flujo másico de la entrada a la columna de destilación
F3 = F3_mass*((aE1_mass/m4) + (mE1_mass/m2) + (hE1_mass/m3)); %kmol/h flujo entrada columna
zF3a = (aE1_mass*F3_mass/m4)/F3;
zF3 = 1 - zF3a; %Fracción molar de m+h

TF3 = T_func(zF3);
HF3 = Fh(zF3);

xD = 1 - aE_mol;
xW = 1 - aB_mol;

%Rmin ---------------------------------------------------------------------

xRmin = zF3;

%Función de la recta de la isoterma
isoT_ex = @(x) interp1([zF3 y_func(xRmin)],[HF3 FH(y_func(xRmin))],x, 'linear','extrap');

%Entalpías de los puntos de diferencias hallados sobre la recta que surge de
%la isoterma
HdeltaD = isoT_ex(xD);
HdeltaW = isoT_ex(xW);

%Cálculo de la reflujo mínimo
Rmin = (HdeltaD - FH(xD))/(FH(xD)-Fh(xD));  

%Valores de x tomados para trazar líneas de operación con Rmin
x_xDmin = linspace(xRmin, xD-0.001, 20);
x_xWmin = linspace(xW+0.001, xRmin, 20);

%Ahora se calcula el y que completa el punto sobre la línea de operación para cada caso
for i = 1 : 20
    Obj2_D = @(y) ((FH(y) - HdeltaD)/(y - xD))-((Fh(x_xDmin(i)) - HdeltaD)/(x_xDmin(i) - xD));
    y_xDmin(i) = fsolve(Obj2_D, 0.5);
    
    Obj2_W = @(y) ((FH(y) - HdeltaW)/(y - xW))-((Fh(x_xWmin(i)) - HdeltaW)/(x_xWmin(i) - xW));
    y_xWmin(i) = fsolve(Obj2_W, 0.1);
end
clc

%Con las parejas de puntos se crean las funciones para las líneas de operación
Op_Wmin = @(x) interp1(x_xWmin, y_xWmin, x);
Op_Dmin = @(x) interp1(x_xDmin, y_xDmin, x);

%% Gráficos destilación

figure('Color','White')
plot(Tb,K(2,Tb),Tb,K(3,Tb),Tb,K(4,Tb))
legend('Metilciclohexano','Heptano','Anilina')
xlabel('Temperaturas, °C')
ylabel(' K value, K = Pvap/Patm')
axis tight
grid minor


%% Gráficos del equilibrio

figure('Color','White')
t = tiledlayout(2,2)  ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile([1,2])
plot(x,Tb,'b',y,Tb,'r')
legend('Temperatura de humbo húmedo','Temperatura de bulbo seco')
xlabel('Fracción molar de metilciclohexano (x_{1}, y_{1})')
ylabel('Temperatura,°C')
xlim([0 1])
grid minor

nexttile 
plot(x,y,'b',[0 1],[0 1],'k')
xlabel('Fracción molar en el líquido, x_{1}')
ylabel('Fracción molar en el vapor ,y_{1}')
grid minor

nexttile
plot(x,HL/1000,'b',y,HV/1000,'r')
hold on
for i = [1:2:16 20:10:100]
    plot([x(i) y_func(x(i))], [Fh(x(i)) FH(y_func(x(i)))]/1000, 'k')
end    
xlim([0 1])
xlabel('Fracción molar de metilciclohexano (x_{1}, y_{1})')
ylabel('Entalpía, kJ/mol')
grid minor

%% Gráficos Rmin destilación

% Linea que pasa por el equlibrio para los valores mínimos
figure('Color','White')
nexttile
hold on
plot(x,HL/1000,'b',y,HV/1000,'r')
for i = [1:2:16 20:10:100]
   plot([x(i) y_func(x(i))], [Fh(x(i)) FH(y_func(x(i)))]/1000,'color', [0.7 0.7 0.7])
end    
xlabel('Fracción molar de mh (x_{1}, y_{1})')
ylabel('Entalpía, kJ/mol')
grid on
plot(zF3, HF3/1000, 'ok','MarkerFaceColor','r')
plot([xW, zF3, xD],isoT_ex([xW, zF3, xD])/1000, 'k')
plot(xW, isoT_ex(xW)/1000, '^k','MarkerFaceColor','b')
plot(xD, isoT_ex(xD)/1000, '^k','MarkerFaceColor','y')

xlim([0 1])


%% Curva de operación mínima
figure('Color','White')
t = tiledlayout(1,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(x,HL/1000,'b',y,HV/1000,'r')
xlabel('Fracción molar de mh (x_{1}, y_{1})')
ylabel('Entalpía, kJ/mol')
grid on
hold on
for i=1:20
    plot([x_xDmin(i) xD],[Fh(x_xDmin(i)) isoT_ex(xD)]/1000, 'color', [0.7 0.7 0.7])
    %plot(y_xD(i),FH(y_xD(i)), 'or')

    plot([y_xWmin(i) xW],[FH(y_xWmin(i)) isoT_ex(xW)]/1000, 'color', [0.5 0.5 0.5])
    %plot(x_xW(i),Fh(x_xW(i)), 'ob')
end
plot(zF3, HF3/1000, 'ok','MarkerFaceColor','r')
plot(xW, isoT_ex(xW)/1000, '^k','MarkerFaceColor','b')
plot(xD, isoT_ex(xD)/1000, '^k','MarkerFaceColor','y')
xlim([0 1])
hold off

nexttile
plot([0 1],[0 1],'b',x,y,'r')
xlabel('Fracción molar en el líquido, x_{1}')
ylabel('Fracción molar en el vapor ,y_{1}')
grid minor
hold on
plot(x_xWmin, Op_Wmin(x_xWmin), 'color', [0.3 0.3 0.3])
plot(x_xDmin, Op_Dmin(x_xDmin), 'color', [0.5 0.5 0.5])
hold off

%% Operatoria con 1.5Rmin --------------------------------------------------------------

R = Rmin*1.5;
%Con R se calcula HdeltaD para este caso
HdeltaD = R*((FH(xD)-Fh(xD))) + FH(xD);

%Ahora se halla la el HdeltaW con el que se cumple la regla de la palanca
Obj11 = @(Hdeltaw) ((HdeltaD - HF3)/(xD-zF3))-((Hdeltaw - HF3)/(xW-zF3));
HdeltaW = fsolve(Obj11, -10);

%También se halla el x en el que se van a encontrar las líneas de operación
Obj12 = @(x) ((HdeltaD - HF3)/(xD-zF3))-((Fh(x) - HdeltaD)/(x - xD));
xcorte = fsolve(Obj12, 0.25);

%Se determinan las parejas de puntos para las nuevas operatorias
%Valores de x tomados para trazar líneas de operación: 
x_xD = linspace(xcorte, xD-0.001, 20);
x_xW = linspace(xW+0.001, xcorte, 20);

%Se calcula el y que completa el punto sobre la línea de operación para
%cada x (vectores x_xD y x_xW)
for i = 1 : 20
    Obj5_D = @(y) ((FH(y) - HdeltaD)/(y - xD))-((Fh(x_xD(i)) - HdeltaD)/(x_xD(i) - xD));
    y_xD(i) = fsolve(Obj5_D, 0.5);
    
    Obj5_W = @(y) ((FH(y) - HdeltaW)/(y - xW))-((Fh(x_xW(i)) - HdeltaW)/(x_xW(i) - xW));
    y_xW(i) = fsolve(Obj5_W, 0.1);
end
clc
%Con las parejas de puntos se crean las funciones para las líneas de
%operación
Op_Wex = @(x) interp1(x_xW, y_xW, x, 'linear','extrap');
Op_Dex = @(x) interp1(x_xD, y_xD, x, 'linear','extrap');
Op_W = @(x) interp1(x_xW, y_xW, x);
Op_D = @(x) interp1(x_xD, y_xD, x);

%Estimación del número de etapas ------------------------------------------

Etx(1) = xD;
Ety(1) = Op_Dex(xD);
Ety(2) = Ety(1);
Etx(2) = x_func(Ety(2));
i = 2;
et = 1;
while Etx(i) > xW
    if Etx(i) > xcorte
        Etx(i+1) = Etx(i);
        Ety(i+1) = Op_D(Etx(i+1));
        Ety(i+2) = Ety(i+1);
        Etx(i+2) = x_func(Ety(i+2));
    elseif Etx > xW+0.001
        Etx(i+1) = Etx(i);
        Ety(i+1) = Op_W(Etx(i+1));
        Ety(i+2) = Ety(i+1);
        Etx(i+2) = x_func(Ety(i+2));
    else
        Etx(i+1) = Etx(i);
        Ety(i+1) = Op_Wex(Etx(i+1));
    end
    i = i + 2;
    et = et + 1;
end    

et = size(Etx,2)/2 ; %Número de etapas (contando rehervidor)

for i=1:1:et
    x_et(i) = Etx(2*i-1); %Para concentración en L0.....
    y_et(i) = Ety(2*i-1); %Para concentración en V1..... 
end    

%% Cálculo de los flujos ----------------------------------------------------

W = F3*(zF3-xD)/(xW-xD) ;
D = F3 - W ;

L(1) = R*D; %L0
V(1) = L(1)+D; %V1

for i=2:(et)
    %Enriquecimiento
    if x_et(i) > xcorte
        L(i) = D * ((xD-y_et(i))/(y_et(i)-x_et(i))) ; %L1 a LF
        V(i) = L(i)+D;    
    %Despojamiento
    else
        L(i) = W * ((y_et(i)-xW)/(y_et(i)-x_et(i))) ;
        V(i) = L(i)-W ;
    end
end

Etapas = 0.5:et-0.5;

figure(10)

nexttile
plot(Etapas,L/1000,'b',Etapas,V/1000,'r')
legend('Liquido','Vapor','Location','northwest')
xlim([0 et])
xlabel('No. Etapa')
ylabel('Flujo molar (kmol/h)')
grid minor

nexttile
plot(Etapas, x_et, 'b', Etapas, y_et, 'r')
legend('Liquido','Vapor','Location','southwest')
xlim([0 et])
xlabel('No. Etapa')
ylabel('Fracción molar de metilciclohexano')
grid minor

%Calor --------------------------------------------------------------------

%Calor condensador
Qc = (HdeltaD - Fh(xD))*D ; %kJ/h
%Calor Rehervidor
Qr = (-HdeltaW + Fh(xW))*W ; %kJ/h

%Cálculo de temperaturas en cada etapa ------------------------------------

Etapas = 0:et;
T_et = T_func(x_et);

figure('Color','White')
plot([Etapas], [T_et T_func(xW)], 'b')
xlim([0 et])
xlabel('No. Etapa')
ylabel('Temperatura (°C)')
grid minor

%% Gráficos con R = 1.5*Rmin ------------------------------------------------

figure('Color','White')
t = tiledlayout(2,1)  ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

plot(x,HL/1000,'b',y,HV/1000,'r')
xlabel('Fracción molar de etanol (x_{1}, y_{1})')
ylabel('Entalpía, kJ/mol')
grid on
hold on
for i=1:20
    plot([x_xD(i) xD],[Fh(x_xD(i)) HdeltaD]/1000, 'color', [0.1 0.3 0.3])
    %plot(y_xD(i),FH(y_xD(i)), 'or')

    plot([y_xW(i) xW],[FH(y_xW(i)) HdeltaW]/1000, 'color', [0.5 0.5 0.5])
    %plot(x_xW(i),Fh(x_xW(i)), 'ob')
end

plot([xW xD], isoT_ex([xW xD])/1000, 'y')

xlim([0 1])
hold off

%% Gráficos de H-x y y-x etapas
figure('Color','White')
t = tiledlayout(2,1)  ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(x,HL/1000,'b',y,HV/1000,'r')
xlabel('Fracción molar de etanol (x_{1}, y_{1})')
ylabel('Entalpía, kJ/mol')
hold on
for i=1:et
    if x_et(i) > xcorte
        plot([xD x_et(i)],[HdeltaD Fh(x_et(i))]/1000,'k')
        plot([x_et(i+1) y_func(x_et(i+1))],[Fh(x_et(i+1)) FH(y_func(x_et(i+1)))]/1000,'g')
    else
        plot([xW y_et(i)],[HdeltaW FH(y_et(i))]/1000, 'k') 
        plot([y_et(i-1) x_func(y_et(i-1))],[FH(y_et(i-1)) Fh(x_func(y_et(i-1)))]/1000,'g')
    end
end

plot([xW xD],[HdeltaW HdeltaD]/1000, 'color',[0.6350 0.0780 0.1840])
plot(zF3, HF3/1000, 'ok','MarkerFaceColor','y')
plot(xW, HdeltaW/1000, '^k','MarkerFaceColor','b')
plot(xD, HdeltaD/1000, '^k','MarkerFaceColor','r')

plot([xW y_et(et)], [Fh(xW) FH(y_et(et))]/1000, 'g')
xlim([0 1])
hold off
grid on

nexttile
plot([0 1],[0 1],'b',x,y,'r',...
    x_xW, Op_W(x_xW), 'k', x_xD, Op_D(x_xD), 'k',Etx, Ety, 'k')
xlabel('Fracción molar en el líquido, x_{1}')
ylabel('Fracción molar en el vapor ,y_{1}')
xlim([0 1])
grid minor

%% y-x, con y sin etapas
figure('Color','white')
t = tiledlayout(2,1) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot([0 1],[0 1],'b',x,y,'r',...
    x_xW, Op_W(x_xW), 'k', x_xD, Op_D(x_xD), 'k')
xlabel('Fracción molar en el líquido, x_{1}')
ylabel('Fracción molar en el vapor ,y_{1}')
xlim([0 1])
grid on

nexttile
plot([0 1],[0 1],'b',x,y,'r',...
    x_xW, Op_W(x_xW), 'k', x_xD, Op_D(x_xD), 'k')
xlabel('Fracción molar en el líquido, x_{1}')
ylabel('Fracción molar en el vapor ,y_{1}')
hold on
plot(Etx, Ety, 'k')
xlim([0 1])
grid on
hold off

%SALIDAS DEL SISTEMA ------------------------------------------------------

S1 = F1*J/1000; %kmol/h
pS1 = y1(1); mS1 = y1(2); hS1 = y1(3);

aC_mass = NC/(1 + NC);
mC_mass = XC*(1-aC_mass);
hC_mass = 1 - aC_mass - mC_mass;
S2_mass = C*(1 + NC); %kg/h Flujo másico en C
S2 = S2_mass*((aC_mass/m4) + (mC_mass/m2) + (hC_mass/m3)); %kmol/h flujo C
aS2 = (aC_mass*S2_mass/m4)/S2;
mS2 = (1 - aS2)*m_h/(1+m_h); %Fracción molar de m en C
hS2 = 1-aS2 - mS2;

aV_mass = NRNp/(1 + NRNp);
mV_mass = XRNp*(1-aV_mass);
hV_mass = 1 - aV_mass - mV_mass;
S3_mass = Ve*(1 + NRNp); %kg/h Flujo másico en V sin propano
S3 = S3_mass*((aV_mass/m4) + (mV_mass/m2) + (hV_mass/m3)) + (F_2*x1(1)/1000); %kmol/h flujo V
aS3 = (Ve*NRNp/m4)/S3;
mS3 = (Ve*XRNp/m2)/S3; %Fracción molar de m en V
hS3 = (Ve*(1-XRNp)/m3)/S3;
pS3 = 1-aS3-mS3-hS3;

%% GRÁFICOS EXTRACCIÓN ----------------------------------------------------

%Gráfico de lineas para construir la operatoria
figure('Color','white')

hold on
plot(X,N_R,'b',Y,N_E,'r')

plot([XDeltaEmin XB XC],[NDeltaEmin NB NC],'g')
plot([XS XRNp XDeltaRmin],[NS NRNp NDeltaRmin],'g')
plot([XDeltaR XF XDeltaE],[NDeltaR NF NDeltaE],'k')
plot([XRNp XDeltaR],[NRNp NDeltaR],'g')

for j = 1:13
    plot([X(j) Y(j)],[N_R(j) N_E(j)],'Color',[0.3 0.3 0.3])
end
for i=1:20
    plot([X_opE(i) Y_opE(i) XDeltaE],[N_R_func(X_opE(i)) N_E_func(Y_opE(i)) NDeltaE], 'color', [0.7 0.7 0.7])
    plot([Y_opR(i) X_opR(i) XDeltaR],[N_E_func(Y_opR(i)) N_R_func(X_opR(i)) NDeltaR], 'color', [0.8 0.8 0.8])
end
grid minor

plot([XDeltaRmin XF XRmin XDeltaEmin],[NDeltaRmin NF N_R_func(XRmin) NDeltaEmin],...
    'color',[0.6350 0.0780 0.1840],'linewidth',1)

legend .
plot(XDeltaE, NDeltaE, 'ok','MarkerFaceColor','r','DisplayName','\DeltaE)')
plot(XDeltaR, NDeltaR, 'ok','MarkerFaceColor','b','DisplayName','\DeltaR')
plot(XDeltaEmin, NDeltaEmin, '^k','MarkerFaceColor','r','DisplayName','\DeltaE_{min}')
plot(XDeltaRmin, NDeltaRmin, '^k','MarkerFaceColor','b','DisplayName','\DeltaR_{min}')
plot(XF, NF, 'ok','MarkerFaceColor','y','DisplayName','F')
plot(XB, NB, 'sk','MarkerFaceColor',[0.9290 0.6940 0.1250],'DisplayName','B')
plot(XC, NC, 'sk','MarkerFaceColor',[0.3010 0.7450 0.9330],'DisplayName','C')
plot(XS, NS, 'sk','MarkerFaceColor',[0.4660 0.6740 0.1880],'DisplayName','S')
xlabel('X,Y (kgM/kgM+H)')
ylabel('N (kgA/kgM+H)')

ylim([-20 50])

%%
figure('Color','white')
t = tiledlayout(2,1) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

%Gráfico de etapas en N-XY
nexttile
plot(X,N_R,'b',Y,N_E,'r')
hold on
for i=1:Et
    if X_et(i) > Xcorte
        plot([XDeltaE X_et(i)],[NDeltaE N_R_func(X_et(i))],'k')
        plot([X_et(i+1) Y_func(X_et(i+1))],[N_R_func(X_et(i+1)) N_E_func(Y_func(X_et(i+1)))],'g')
    else
        plot([XDeltaR Y_et(i)],[NDeltaR N_E_func(Y_et(i))], 'k') 
        plot([Y_et(i-1) X_func(Y_et(i-1))],[N_E_func(Y_et(i-1)) N_R_func(X_func(Y_et(i-1)))],'g')
    end
end
plot([XDeltaR XRNp], [NDeltaR N_R_func(XRNp)], 'k')
plot([XRNp Y_et(Et)], [N_R_func(XRNp) N_E_func(Y_et(Et))], 'g')
plot(XDeltaE, NDeltaE, '^k','MarkerFaceColor','r','DisplayName','\DeltaE)')
plot(XDeltaR, NDeltaR, '^k','MarkerFaceColor','b','DisplayName','\DeltaR')
plot(XF, NF, 'ok','MarkerFaceColor','y','DisplayName','F')
xlabel('X,Y (kgM/kgM+H)')
ylabel('N (kgA/kgM+H)')

grid minor
ylim([-25 35])

nexttile 
plot([0 1],[0 1],'b',X_opR, Op_R(X_opR), 'k', X_opE, Op_E(X_opE), 'k', X,Y,'r',EtX, EtY, 'k')
xlabel('X(kgM/kgM+H)')
ylabel('X(kgM/kgM+H)')
grid minor

%% Gráfico de etapas en X-Y
figure('Color','white')
t = tiledlayout(1,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot([1:Et+1], E_et, 'r', [1:Et+1], R_et, 'b')
legend('Extracto','Refinado')
grid minor
xlabel('No. Etapa')
ylabel('Flujo másico (Libre de solvente), kg')

nexttile
plot([1:Et+1], Y_et, 'r', [1:Et+1], X_et, 'b')
grid minor
legend('Extracto','Refinado')
xlabel('No. Etapa')
ylabel('X,Y (kgM/kgM+H)')

%FUNCIONES ----------------------------------------------------------------

%Coeficiente de distribución - con Antoine
function f = K(i, T)
    global P;
    A = [6.80398 6.823 6.89385 7.2418];
    B = [803.81 1270.763 1264.37 1675.3];
    C = [246.99 221.416 216.636 200.01];
    
    f = (10.^(A(i) - B(i)./(T + C(i))))*101.325/(760*P);
end
