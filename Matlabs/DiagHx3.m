clear all
clc

%Datos del ejercicio
R = 8.3145 ; %J/molK
P = 760 ;%mmHg
Tc1 = 513.92 ; %K
Tc2 = 647.14 ; %K
T_ref = 25 ; %°C
M1 = 46 ; %g/mol
M2 = 18 ; %g/mol

%Coeficiente de actividad (Van Laar dos parametros)
A12 = 1.6798 ; A21 = 0.9227 ;

%Equilibrio de Fases
A1 = 7.68117 ; B1 = 1332.04 ; C1 = 199.200 ;
A2 = 8.07131 ; B2 = 1730.63 ; C2 = 233.426 ;

%Presión de equilibrio
x1 = linspace(0,1,51) ;
x1 = x1' ;
[gama1, gama2] = VanLaar(A12,A21,x1) ;

for i=1:size(x1,1)
f = @(t) x1(i)*gama1(i)*Pvap(A1,B1,C1,t) + (1-x1(i))*gama2(i)*Pvap(A2,B2,C2,t)-P ;
T(i,1) = fsolve(f,50) ;
end

y1 = x1.*gama1.*Pvap(A1,B1,C1,T)./P ;

%Composición y temperatura del azeotrópo
x_AZ = fsolve(@(x) interp1(x1,y1,x)-x, 0.9) ;
T_AZ = interp1(x1,T,x_AZ) ;
clc

figure('Color','White')
t=tiledlayout(2,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot([0 1],[0 1],'b',x1,y1,'r',x_AZ,x_AZ,'oy')
xlabel('Liquid Mole Fraction, x_{1}')
ylabel('Vapor Mole Fraction ,y_{1}')
grid minor

nexttile([2 1])
plot(x1,T,'b',y1,T,'r',x_AZ,T_AZ,'oy')
xlabel('Ethanol Mole Fraction (x_{1}, y_{1})')
ylabel('Temperature, °C')
grid minor

Tx = linspace(T(1),T(end),50)' ;
[Tx, ind] = sort([Tx ;T_AZ]) ;
[v indAZ] = max(ind) ;
Ty = Tx ;

x1 = interp1(T,x1,Tx) ;
[gama1, gama2] = VanLaar(A12,A21,x1) ;
y1 = x1.*gama1.*Pvap(A1,B1,C1,Ty)./P ;

%Capacidades calorificas:
CpL2 = 75.66 ; CpL1 = 163.8 ;
CpV2 = 38.29 ; CpV1 = 73.98 ;

%Entalpia en exceso:
T_K = Tx + 273.15 ; %K
a0 =-3.63868e5 + 1.83829e3.*T_K - 2.32763.*T_K.^2 ; 
a1 = 9.25982e5 - 4.83586e3.*T_K + 6.37228.*T_K.^2 ;
a2 = -14.04894e5 + 7.51661e3.*T_K - 10.11280.*T_K.^2 ;
a3 = 10.91318e5 - 5.89498e3.*T_K + 7.98868.*T_K.^2 ;
a4 = -2.79986e5 + 1.50557e3.*T_K - 2.03127.*T_K.^2  ;
HE = x1.*(1-x1).*(a0 + a1.*x1.^0.5+a2.*x1.^1.5+a3.*x1.^2.5+a4.*x1.^4.5) ; %J/mol

%Calores de vaporización:
T_refK = T_ref+273.15 ;
C11 = 5.6900e7 ; C12 = 0.3359 ; C13 = 0 ; C14 = 0 ;
C21 = 5.2053e7 ; C22 = 0.3199 ; C23 =-0.212 ; C24 = 0.25795;
lambda1 = lambda(C11,C12,C13,C14,T_refK,Tc1) ;
lambda2 = lambda(C21,C22,C23,C24,T_refK,Tc2) ;

%Entalpias de equilibrio
HL = x1.*CpL1.*(Tx-T_ref) + (1-x1).*CpL2.*(Tx-T_ref) + HE ; %kJ/kmol
HV = y1.*(lambda1 + CpV1.*(Ty-T_ref)) + (1-y1).*(lambda2 + CpV2.*(Ty-T_ref)) ; %kJ/kmol

%Masa molecular promedio
Mpromx = x1.*M1 + (1-x1).*M2 ;
Mpromy = y1.*M1 + (1-y1).*M2 ;

Hl = HL./Mpromx/4.18; %kcal/kg
Hv = HV./Mpromy/4.18 ; %kcal/kg

wx1 = x1./Mpromx*M1 ;
wy1 = y1./Mpromy*M1 ;

nexttile
x = [x1 y1] ;
Hm = [HL HV]/1000 ;
H = [Hl Hv]*4.18 ; %kJ/kg
plot(x1,HL/1000,'b',y1,HV/1000,'r')
xlabel('Fracción molar de Etanol (x_{1}, y_{1})')
ylabel('Entalpia, kJ/mol')
axis([0 1 0 60])
grid on
hold on
for i=[1:2:10 11:4:size(x1,1)-4 50]
  plot(x(i,:),Hm(i,:),'k')
end
 plot(x(indAZ,:),H(indAZ,:)/1000,'y')

figure('Color','White')
plot(wx1,Hl,'b',wy1,Hv,'r')
axis([0 1 0 700])
xlabel('Ethanol Weight Fraction')
ylabel('Enthalpy, kCal/kg')
grid on

%Isotermas de equilibrio
xw = [wx1 wy1] ;
Hw = [Hl Hv] ;
hold on
for i= [1:2:10 11:4:size(x1,1)-4 50]
  plot(xw(i,:),Hw(i,:),'k')
  txt = sprintf(' %.2f °C ',Tx(i)) ;
  posx = mean(xw(i,:),2) ;
  posy = mean(Hw(i,:),2) ;
  angle = atan((Hw(i,2)-Hw(i,1))/(abs(xw(i,2)-xw(i,1))*1000))*180/pi ;
  text(posx,posy,txt,'HorizontalAlignment','center',...
      'Rotation',angle,'BackgroundColor','w')
end

%Dew linew
  posy = mean([Hw(1,2) Hw(end,2)],2) ;
  angle = atan((Hw(1,2)-Hw(end,2))/((xw(1,1)-xw(end,1))*1000))*180/pi ;
  text(0.5,posy,'Dew Line','HorizontalAlignment','center',...
      'Rotation',angle,'BackgroundColor','w','Color','r')
 %Boiling line
 posy = mean([Hw(1,1) Hw(end,1)],2) ;
  angle = atan((Hw(1,1)-Hw(end,1))/((xw(1,1)-xw(end,1))*1000))*180/pi ;
  text(0.5,posy,'Boiling Line','HorizontalAlignment','center',...
      'Rotation',angle,'BackgroundColor','w','Color','b')
  
  
function [gama1,gama2]=VanLaar(A12,A21,xA)
gama1 = exp(A12.*(A21.*(1-xA)./(A12.*xA + A21.*(1 - xA))).^2) ;
gama2 = exp(A21.*(A12.*xA./(A12.*xA + A21.*(1-xA))).^2) ;
end
  
function f=lambda(C1,C2,C3,C4,T,Tc) 
    f =10^-3*C1.*(1-T/Tc).^(C2+C3.*T/Tc+C4.*(T/Tc).^2) ;
end

  function f = Pvap(A,B,C,T)
    f = 10.^(A - B./(T + C)) ;
end