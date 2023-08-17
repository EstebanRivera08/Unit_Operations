clear all
clc

%1-Etanol
%2-Agua

%Datos del ejercicio
Rgases = 8.3145 ; %J/molK
P = 760 ;%mmHg
Tc1 = 513.92 ; %K
Tc2 = 647.14 ; %K
T_ref = 0 ; %°C
M1 = 46 ; %g/mol
M2 = 18 ; %g/mol

%Coeficiente de actividad (Van Laar dos parametros)
A12 = 1.6798 ; A21 = 0.9227 ;

%Equilibrio de Fases - Constantes ec de Antoine
A1 = 7.68117 ; B1 = 1332.04 ; C1 = 199.200 ;
A2 = 8.07131 ; B2 = 1730.63 ; C2 = 233.426 ;

%Presión de equilibrio 
x1 = linspace(0,1,101) ;  %Valores de composición en líquido con los que se va a trabajar
x1 = x1' ;
[gama1, gama2] = VanLaar(A12,A21,x1) ; %Coeficientes de actividad en cada caso

%Determinación de la temperatura para cada caso con la que se cumple la
%presión con la ley de Raoult modificada (T de equilibrio)
for i=1:size(x1,1) 
f = @(t) x1(i)*gama1(i)*Pvap(A1,B1,C1,t) + (1-x1(i))*gama2(i)*Pvap(A2,B2,C2,t)-P ; %Solver para T
Tx(i,1) = fsolve(f,50) ; %T del equilibrio para cada valor de xi
end

%A partir de la deficinición de equilibrio con la fugacidad y haciendo
%algunas modificaciones, se calcula cada composicón en el vapor en
%equilibrio con las dadas para xi
y1 = x1.*gama1.*Pvap(A1,B1,C1,Tx)./P ; %Vector
%Y se disponen estos valores como una función:
y_func = @(x) interp1(x1, y1, x); %Función
x_func = @(y) interp1(y1, x1, y); %Función para x en términos de y 
T_func = @(x) interp1(x1, Tx, x); %Función de T para cada x en equilibrio

%Composición y temperatura del azeótropo
x_AZ = fsolve(@(x) interp1(x1,y1,x)-x, 0.9) ;
T_AZ = interp1(x1,Tx,x_AZ) ;
clc

%Capacidades calorificas:
CpL2 = 75.66 ; CpL1 = 163.8 ;
CpV2 = 38.29 ; CpV1 = 73.98 ;

%Calores de vaporización:
T_refK = T_ref+273.15 ;
C11 = 5.6900e7 ; C12 = 0.3359 ; C13 = 0 ; C14 = 0 ;
C21 = 5.2053e7 ; C22 = 0.3199 ; C23 =-0.212 ; C24 = 0.25795;
lambda1 = lambda(C11,C12,C13,C14,T_refK,Tc1) ;
lambda2 = lambda(C21,C22,C23,C24,T_refK,Tc2) ;

%Entalpias de equilibrio L-líquido, V-vapor
HL = (x1.*CpL1.*(Tx-T_ref) + (1-x1).*CpL2.*(Tx-T_ref) + HE(Tx, x1))/1000 ; %kJ/mol
HV = (y1.*(lambda1 + CpV1.*(Tx-T_ref)) + (1-y1).*(lambda2 + CpV2.*(Tx-T_ref)))/1000; %kJ/mol

Fh = @(x) interp1(x1,HL,x);
FH = @(x) interp1(y1,HV,x);

%--------------------------------------------------------------------------

%VARIABLES DADAS

zF = 0.2;
xD = 0.76;
xW = 0.02;
F = 1000*453.59; %mol/h
TF = 75; %°C
W = F*(zF-xD)/(xW-xD) ;
D = F - W ;

%Función para la entalpía en líquido comprimido
hL = @(x, Temp) (x*CpL1*(Temp-T_ref) + (1-x)*CpL2*(Temp-T_ref)+ HE(Temp, x))/1000 ; %kJ/mol

hF = hL(zF, TF); %Entalpía de la corriente de entrada

%--------------------------------------------------------------------------

%PRIMERA APROXIMACIÓN Rmin
%Se extiende una recta isoterma de la región liq-vap hasta tocar a la
%corriente entrada. Se busca el caso en el que la concentración del liq-sat
%tenga esta temperatura (Se halla resolviendo la función objetivo)
Obj = @(x) ((FH(y_func(x)) - Fh(x))/(y_func(x) - x)) - ((hF - Fh(x))/(zF - x));
xRmin = fsolve(Obj, 0.5);
%Función de la recta de la isoterma
isoT = @(x) interp1([zF y_func(xRmin)],[hF FH(y_func(xRmin))],x);
isoT_ex = @(x) interp1([zF y_func(xRmin)],[hF FH(y_func(xRmin))],x, 'linear','extrap');

%Entalpías de los puntos de diferencias hallados sobre la recta que surge de
%la isoterma
HdeltaD = isoT_ex(xD);
HdeltaW = isoT_ex(xW);

HdeltaWmin = HdeltaW ;
HdeltaDmin = HdeltaD ;

%Cálculo de la reflujo mínimo
Rmin = (HdeltaD - FH(xD))/(FH(xD)-Fh(xD));

%CÁLCULO DE Rmin
%Valores de x tomados para trazar líneas de operación 
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

%Con las parejas de puntos se crean las funciones para las líneas de
%operación
Op_Wmin = @(x) interp1(x_xWmin, y_xWmin, x);
Op_Dmin = @(x) interp1(x_xDmin, y_xDmin, x);

%Gráficos Rmin
t=tiledlayout(1,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
x = [x1 y1] ;
Hm = [HL HV]; %kJ/mol
plot(x1,HL,'b',y1,HV,'r')
xlabel('Fracción molar de etanol (x_{1}, y_{1})')
ylabel('Entalpía, kJ/mol')
grid minor
hold on
plot([xW xD],[HdeltaW HdeltaD], 'k','LineWidth',1.1)
for i=[1:4 4:2:16 16:8:80 80:4:100]
  plot(x(i,:),Hm(i,:),'k')
end

plot(zF, hF, 'ok','MarkerFaceColor','m')
plot(xRmin,Fh(xRmin),'ok','MarkerFaceColor','b')
plot(y_func(xRmin),FH(y_func(xRmin)),'ok','MarkerFaceColor','r')
plot(xW,HdeltaW,'^k','MarkerFaceColor','y')
plot(xD,HdeltaD,'^k','MarkerFaceColor','y')
hold off

nexttile
plot([0 1],[0 1],'b',x1,y1,'r',x_AZ,x_AZ,'oy')
xlabel('Fracción molar en el líquido, x_{1}')
ylabel('Fracción molar en el vapor ,y_{1}')
grid minor
hold on
plot(x_xWmin, Op_Wmin(x_xWmin), 'k')
plot(x_xDmin, Op_Dmin(x_xDmin), 'k')
hold off

%% --------------------------------------------------------------------------

%CON 1.5Rmin

R = Rmin*1.5;
%Con R se calcula HdeltaD para este caso
HdeltaD = R*((FH(xD)-Fh(xD))) + FH(xD);

%Ahora se halla la el HdeltaW con el que se cumple la regla de la palanca
Obj3 = @(Hdeltaw) ((HdeltaD - hF)/(xD-zF))-((Hdeltaw - hF)/(xW-zF));
HdeltaW = fsolve(Obj3, -10);

%También se halla el x en el que se van a encontrar las líneas de operación
Obj4 = @(x) ((HdeltaD - hF)/(xD-zF))-((Fh(x) - HdeltaD)/(x - xD));
xcorte = fsolve(Obj4, 0.25);

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

%Estimación del número de etapas
EtX(1) = xD;
EtY(1) = Op_Dex(xD);
EtY(2) = EtY(1);
EtX(2) = x_func(EtY(2));
i = 2;
Et = 1;
while EtX(i) > xW
    if EtX(i) > xcorte
        EtX(i+1) = EtX(i);
        EtY(i+1) = Op_D(EtX(i+1));
        EtY(i+2) = EtY(i+1);
        EtX(i+2) = x_func(EtY(i+2));
    elseif EtX > xW+0.001
        EtX(i+1) = EtX(i);
        EtY(i+1) = Op_W(EtX(i+1));
        EtY(i+2) = EtY(i+1);
        EtX(i+2) = x_func(EtY(i+2));
    else
        EtX(i+1) = EtX(i);
        EtY(i+1) = Op_Wex(EtX(i+1));
    end
    i = i + 2;
    Et = Et + 1;
end    

Et = size(EtX,2)/2 ; %Número de etapas (contando rehervidor)

for i=1:1:Et
    x_et(i) = EtX(2*i-1); %Para concentración en L0.....L10
    y_et(i) = EtY(2*i-1); %Para concentración en V1.....V11  
end    
%
%Gráficos con R = 1.5*Rmin
%
% Gráficos de H-x y y-x etapas

figure()
t = tiledlayout(2,1) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(x1,HL,'b',y1,HV,'r',xD,HdeltaD,'oy',xW,HdeltaW,'oy')
xlabel('Fracción molar de etanol (x_{1}, y_{1})')
ylabel('Entalpía, kJ/mol')
hold on
for i=1:Et
    if x_et(i) > xcorte
        plot([xD x_et(i)],[HdeltaD Fh(x_et(i))],'k')
        plot([x_et(i+1) y_func(x_et(i+1))],[Fh(x_et(i+1)) FH(y_func(x_et(i+1)))],'g')
    else
        plot([xW y_et(i)],[HdeltaW FH(y_et(i))], 'k') 
        plot([y_et(i-1) x_func(y_et(i-1))],[FH(y_et(i-1)) Fh(x_func(y_et(i-1)))],'g')
    end
end
plot([xW xW], [HdeltaW Fh(xW)], 'k',[xW y_et(Et)], [Fh(xW) FH(y_et(Et))], 'g')
plot([xW xD],[HdeltaW HdeltaD],'Color',[0.5 0.3 0.2])
plot(zF, hF, 'ok','MarkerFaceColor','y')
plot(xW,HdeltaW,'ok','MarkerFaceColor','b')
plot(xD,HdeltaD,'ok','MarkerFaceColor','r')

hold off
grid minor

nexttile
plot([0 1],[0 1],'b',x1,y1,'r',x_AZ,x_AZ,'oy',...
    x_xW, Op_W(x_xW), 'k', x_xD, Op_D(x_xD), 'k',EtX, EtY, 'k')
xlabel('Fracción molar en el líquido, x_{1}')
ylabel('Fracción molar en el vapor ,y_{1}')
grid minor

% y-x, con y sin etapas
figure()
t=tiledlayout(1,2)
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
plot([0 1],[0 1],'b',x1,y1,'r',x_AZ,x_AZ,'oy',...
    x_xW, Op_W(x_xW), 'k', x_xD, Op_D(x_xD), 'k')
xlabel('Fracción molar en el líquido, x_{1}')
ylabel('Fracción molar en el vapor ,y_{1}')
grid on

nexttile
plot([0 1],[0 1],'b',x1,y1,'r',x_AZ,x_AZ,'oy',...
    x_xW, Op_W(x_xW), 'k', x_xD, Op_D(x_xD), 'k')
xlabel('Fracción molar en el líquido, x_{1}')
ylabel('Fracción molar en el vapor ,y_{1}')
hold on
plot(EtX, EtY, 'k')
grid on
hold off

%--------------------------------------------------------------------------
%Cálculo de los flujos

L(1) = R*D; %L0
V(1) = L(1)+D; %V1

for i=2:(Et)
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

Etapas = 0.5:Et-0.5;

figure()

nexttile
plot(Etapas,L/1000,'b',Etapas,V/1000,'r')
legend('Liquido','Vapor','Location','northwest')
xlim([0 Et])
xlabel('No. Etapa')
ylabel('Flujo molar (kmol/h)')
grid minor

nexttile
plot(Etapas, x_et, 'b', Etapas, y_et, 'r')
legend('Liquido','Vapor','Location','southwest')
xlim([0 Et])
xlabel('No. Etapa')
ylabel('Fracción molar de Etanol')
grid minor

%Calor condensador
Qc = (HdeltaD - Fh(xD))*D ; %kJ/h
%Calor Rehervidor
Qr = (-HdeltaW + Fh(xW))*W ; %kJ/h

%--------------------------------------------------------------------------
%Cálculo de temperaturas en cada etapa
Etapas = 0:Et;
T_et = T_func(x_et);

figure()
plot([Etapas], [T_et T_func(xW)], 'b')
xlim([0 Et])
xlabel('No. Etapa')
ylabel('Temperatura (°C)')
grid minor
%
%--------------------------------------------------------------------------
%Resultados

sprintf('Qr: %.2f, kJ/h', Qr/3600)
sprintf('Qc: %.2f, kJ/h', Qc/3600)
sprintf('Num etapas: %.2f', Et)
sprintf('W: %.2f, kmol/h', W/1000)
sprintf('D: %.2f, kmol/h', D/1000)
sprintf('L1: %.2f, kmol/h', L(2)/1000)
sprintf('V2: %.2f, kmol/h', V(2)/1000)
sprintf('Ln-1: %.2f, kmol/h', L(Et)/1000)
sprintf('Vn: %.2f, kmol/h', V(Et)/1000)


%--------------------------------------------------------------------------
%Flujos de vapor de calentamiento y agua de enfriamiento

%Vap agua entra a 180 y sale a 110 °C
Cp_vap = 2.080; %kJ/kg*K
m_vap = Qr/(Cp_vap*(180-110)*3600); %kg/s

%Liq agua entra a 35 y sale a 43 °C
Cp_liq = 4.1813; %kJ/kg*K
m_liq = Qc/(Cp_liq*(43-35)*3600); %kg/s

sprintf('M vap: %.2f, kg/s', m_vap)
sprintf('M liq: %.2f, kg/s', m_liq)


%%


function exceso = HE(Tx, x1) %Función para determinar la entalpía en exceso
    T_K = Tx + 273.15 ; %K
    a0 =-3.63868e5 + 1.83829e3.*T_K - 2.32763.*T_K.^2 ; 
    a1 = 9.25982e5 - 4.83586e3.*T_K + 6.37228.*T_K.^2 ;
    a2 = -14.04894e5 + 7.51661e3.*T_K - 10.11280.*T_K.^2 ;
    a3 = 10.91318e5 - 5.89498e3.*T_K + 7.98868.*T_K.^2 ;
    a4 = -2.79986e5 + 1.50557e3.*T_K - 2.03127.*T_K.^2  ;
    exceso = x1.*(1-x1).*(a0 + a1.*x1.^0.5+a2.*x1.^1.5+a3.*x1.^2.5+a4.*x1.^4.5) ; %J/mol
end

function [gama1,gama2]=VanLaar(A12,A21,xA) %Coeficientes de actividad
gama1 = exp(A12.*(A21.*(1-xA)./(A12.*xA + A21.*(1 - xA))).^2) ;
gama2 = exp(A21.*(A12.*xA./(A12.*xA + A21.*(1-xA))).^2) ;
end
  
function f=lambda(C1,C2,C3,C4,T,Tc) %Constantes para calor de vaporización 
    f =10^-3*C1.*(1-T/Tc).^(C2+C3.*T/Tc+C4.*(T/Tc).^2) ;
end

  function f = Pvap(A,B,C,T)  %Función Antoine
    f = 10.^(A - B./(T + C)) ;
end