
clear
clc

% Cálculo general

%Dimensiones

global Stoi_Coef R dt

dt = 0.15 ; % m
zF = 15 ; %m
Tf = 1200 ; %K
G = 54.35 ; % kg/m2s PARA ETANO CON IMPUREZAS
%G = 53.52 ; % kg/m2s PARA ETANO PURO
Ntube = 6 ;

% Kinetics

%constants
R = 8.3145 ; %kJ/kmolK

%Global Rates
Stoi_Coef = [ 0 , 0 , 1 ,-1 , 0 , 0 , 0 , 1 , 0 ;...
              0 , 0 ,-1 , 1 , 0 , 0 , 0 ,-1 , 0 ;...
              1 , 0 , 0 ,-2 , 0 , 1 , 0 , 0 , 0 ;...
              1 , 1 , 0 , 0 ,-1 , 0 , 0 , 0 , 0 ;...
             -1 ,-1 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ;...
              0 ,-1 ,-1 , 0 , 0 , 0 , 1 , 0 , 0 ;...
              1 , 0 ,-1 ,-1 , 1 , 0 , 0 , 0 , 0 ] ;
      

% [  CH4    C2H2     C2H4   C2H6    C3H6     C3H8    C4H8     H2      H2O]
% [Metano Acetileno Etileno Etano Propileno Propano Buteno Hidrogeno Agua]
% [   1       2        3      4       5        6      7        8       9 ]

% Condiciones iniciales
A = pi*dt^2/4 ;
m3 = 69*1000 ; %ton/año
P0 = 2.99*101.325 ; %kPa

T0 = 680 + 273.15 ; %°K

xmax = 0.60 ; %Conversión del etano
PF = 1.2*101.325 ; %atm

y3 = 0.01 ; %mol/mol
y4 = 0.982 ; 
y5 = 0.008 ;
mH2O_mEt = 0.4 ; %kgH2O/kgEtano
nH2O_nEt = mH2O_mEt*30.07/18 ; %molH2O/molEtano

%IMPUREZAS
F0 = [0, 0, y3, y4, y5, 0, 0, 0, nH2O_nEt*y4]*G*A/30.07 ; % kmol/s

%SIN IMPUREZAS
%F0 = [0, 0, 0, 1, 0, 0, 0, 0, nH2O_nEt]*G*A/30.07 ; % kmol/s

Mm0 = G*A/sum(F0) ;
l_Mm0 = 1/Mm0 ; %kmol/s
 
iter = 0 ;
dif = 10 ;


while dif > 1e-3 || iter < 20
%Balance de energía
% [  CH4    C2H2     C2H4   C2H6    C3H6     C3H8    C4H8     H2      H2O]
% [Metano Acetileno Etileno Etano Propileno Propano Buteno Hidrogeno Agua]
% [   1       2        3      4       5        6      7        8       9 ]

Cpi = Cp((T0+Tf)/2) ; % kJ/kmolK

sigma = 2.041*10^(-7) ; %Constante de boltzman kJ/h.m^2.K^4
emisiv = 0.835; %Relative effectiveness factor of the bank
F = 0.97 ; %Factor de intercambio
hconv = 30.66; % Film convective heat transfer coefficient; (kJ/h.m2.c)

Tg1 = Tg(Tf,zF,Ntube); %K

q = @(z,T) (hconv*(Tg1-T) + sigma*emisiv*F*(Tg1^4-(T).^4))/60/60 ; %Función de calor kJ/sm2


Tref = 25+273.15 ; % K

Hf0 = [-74.8 , 227.1 , 52.53 , -84.2 , 20.41 , -104.7, 110 , 0 , -241.83]*1000 ; %kJ/kmol

Hf = @(T) Hf0 + Cpi.*(T - Tref); %Entalpías de reacción 1X9 %kJ/kmol

Hr = @(T) Stoi_Coef*Hf(T)'  ; %Vector 7x1 kJ/kmol

dTdz = @(z,T,P,Fi) 1./(Fi*Cpi')*(q(z,T)*pi*dt - A*(ri(T,P,Fi)*Hr(T))) ;

%Balance de energía mecánica (ERGUN)
alfa = 1 ; %m/s , Pa

mu = @(T,Fi) Visc(T,Fi) ; % Pa*s = kg/sm
Re = @(T,Fi) dt*G./mu(T,Fi) ; 
f =  @(T,Fi) 0.046*Re(T,Fi)^(-0.2) ;


d1_Mmdz = @(z,Fi,T,P) 1/(G*A)*(dFidz(1,z,T,P,Fi)+dFidz(2,z,T,P,Fi)+...
    dFidz(3,z,T,P,Fi)+dFidz(4,z,T,P,Fi)+dFidz(5,z,T,P,Fi)+...
    dFidz(6,z,T,P,Fi)+dFidz(7,z,T,P,Fi)+dFidz(8,z,T,P,Fi)+dFidz(9,z,T,P,Fi)); %kmol/kgm

dPdz = @(z,T,P,Fi,l_Mm)  (d1_Mmdz(z,Fi,T,P) + l_Mm.*(1./T.*dTdz(z,T,P,Fi)+...
    2*f(T,Fi)/dt))./(l_Mm./P - (P./(G^2.*R.*T))*1000) ; %kPa/m 

SistemaEc = @(z,T,P,F1,F2,F3,F4,F5,F6,F7,F8,F9,l_Mm) [dFidz(1,z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9]);...
    dFidz(2,z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9]);...
    dFidz(3,z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9]);...
    dFidz(4,z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9]);...
    dFidz(5,z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9]);...
    dFidz(6,z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9]);...
    dFidz(7,z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9]);...
    dFidz(8,z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9]);...
    dFidz(9,z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9]);...
    dTdz(z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9]);...
    d1_Mmdz(z,[F1,F2,F3,F4,F5,F6,F7,F8,F9],T,P);...
    dPdz(z,T,P,[F1,F2,F3,F4,F5,F6,F7,F8,F9],l_Mm)] ;

n = 500 ;
h = 60/n ;
z = linspace(0,60,n+1)' ;
CondIniciales = [F0'; T0; l_Mm0; P0] ;
X(1,:) = CondIniciales' ;
Y(1,:) = X(1,:) ;

for i = 1:size(z,1)-1
   k1 = SistemaEc(z(i),X(i,10),X(i,12),X(i,1),X(i,2),X(i,3),X(i,4),...
        X(i,5),X(i,6),X(i,7),X(i,8),X(i,9),X(i,11))' ;
    
   k2 = SistemaEc(z(i)+h/2,X(i,10)+h/2*k1(10),X(i,12)+h/2*k1(12),...
       X(i,1)+h/2*k1(1),X(i,2)+h/2*k1(2),X(i,3)+h/2*k1(3),X(i,4)+h/2*k1(4),...
       X(i,5)+h/2*k1(5),X(i,6)+h/2*k1(6),X(i,7)+h/2*k1(7),X(i,8)+h/2*k1(8),...
       X(i,9)+h/2*k1(9),X(i,11)+h/2*k1(11))' ;

   k3 = SistemaEc(z(i)+h/2,X(i,10)+h/2*k2(10),X(i,12)+h/2*k2(12),...
       X(i,1)+h/2*k2(1),X(i,2)+h/2*k2(2),X(i,3)+h/2*k2(3),X(i,4)+h/2*k2(4),...
       X(i,5)+h/2*k2(5),X(i,6)+h/2*k2(6),X(i,7)+h/2*k2(7),X(i,8)+h/2*k2(8),...
       X(i,9)+h/2*k2(9),X(i,11)+h/2*k2(11))' ;   

   k4 = SistemaEc(z(i)+h,X(i,10)+h*k3(10),X(i,12)+h*k3(12),...
       X(i,1)+h*k3(1),X(i,2)+h*k3(2),X(i,3)+h*k3(3),X(i,4)+h*k3(4),...
       X(i,5)+h*k3(5),X(i,6)+h*k3(6),X(i,7)+h*k3(7),X(i,8)+h*k3(8),...
       X(i,9)+h*k3(9),X(i,11)+h*k3(11))' ; 
   
   X(i+1,:) = X(i,:) + h/6*(k1+2*k2+2*k3+k4) ;
   
%    X(i+1,:) = X(i,:) + h*SistemaEc(z(i),X(i,10),X(i,12),X(i,1),X(i,2),X(i,3),X(i,4),...
%         X(i,5),X(i,6),X(i,7),X(i,8),X(i,9),X(i,11))' ;
end

%Resultados Euler
% Fi_E = X(:,1:9)*1000 ; %mol
% T_E = X(:,10) ;
% P_E = X(:,12) ;
%x_E = (Fi_E(1,4)-Fi_E(:,4))/Fi_E(1,4) ;

%Resultados Runge-Kutta
Fi = X(:,1:9)*1000 ; %mol
T = X(:,10) ;
P = X(:,12) ;
x= (Fi(1,4)-Fi(:,4))/Fi(1,4) ;

%Funciones
xA = @(y) interp1(z,x,y) ;
zF = fsolve(@(y) xA(y)-xmax,5) ;
clc
F_Etileno = @(y) interp1(z,Fi(:,3),y) ;
Tope = @(y) interp1(z,T,y) ;


dif = abs((Tf-Tope(zF))) ;
Tf = Tope(zF) ;
F_anual = F_Etileno(zF)*28.05/1000/1000*60*60*24*365 ; % Flujo etileno ton/año

Ntube = round(m3/F_anual) ;
F_total = F_anual*Ntube ;

iter = iter +1 ;
end

%% Selectividad y Rendimiento

% [  CH4    C2H2     C2H4   C2H6    C3H6      C3H8     C4H6      H2      H2O]
% [Metano Acetileno Etileno Etano Propileno Propano Butadieno Hidrogeno Agua]
% [   1       2        3      4       5        6        7        8       9 ]

F1 = @(y) interp1(z,Fi(:,1),y) ;
F2 = @(y) interp1(z,Fi(:,2),y) ;
F3 = @(y) interp1(z,Fi(:,3),y) ;
F4 = @(y) interp1(z,Fi(:,4),y) ;
F5 = @(y) interp1(z,Fi(:,5),y) ;
F6 = @(y) interp1(z,Fi(:,6),y) ;
F7 = @(y) interp1(z,Fi(:,7),y) ;
F8 = @(y) interp1(z,Fi(:,8),y) ;
F9 = @(y) interp1(z,Fi(:,9),y) ;

L = linspace(0,zF,100) ;

Rend = F3(L)./(F4(0)-F4(L)) ;
Select = F3(L)./(F1(L)+F2(L)+F5(L)+F6(L)+F7(L)+F8(L)) ;

figure('Color','White')
t =tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(L,Rend,'r')
xlabel('z, m')
ylabel('Y_{Etileno}')
ylim([0.8 1])
grid minor

[x,y] = mesh(X(1:10,1),X(1:10,2)) ;

nexttile
plot(L,Select,'b')
xlabel('z, m')
ylabel('S_{Etileno}')
grid minor

%% Flujos molares
figure('Color','White')
hold on
for j=1:9
plot(z,Fi(:,j))
end
ylim([0 0.75])
xlim([0 zF])
grid minor
ylabel('F_i , mol/s')
xlabel('z, m')
legend('Metano', 'Acetileno', 'Etileno', 'Etano', 'Propileno', 'Propano', 'Butadieno', 'Hidrogeno', 'Agua')

%% Conversión-Temperatura-Presión

figure('Color','White')
t =tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
 
nexttile([2;1])
plot(z,x,'k')
xlim([0 zF])
grid minor
ylabel('x_{Etano}')
xlabel('z, m')

nexttile
plot(z,T,'r')
xlim([0 zF])
grid minor
ylabel('T , K')
xlabel('z, m')

nexttile
plot(z,P,'b')
xlim([0 zF])
grid minor
ylabel('P , kPa')
xlabel('z, m')

%% Masa molecular - Densidad - Viscosidad

Mm = 1./X(:,11) ;
rho = Mm.*P./(R.*T) ;
for i=1:size(Fi,1) 
mu1(i) = Visc(T(i),Fi(i,:)) ;
end

figure('Color','White')
t =tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile([2;1])
plot(z,rho,'k')
xlim([0 zF])
grid minor
ylabel('\rho , kg/m^3')
xlabel('z, m')

nexttile
plot(z,Mm,'b')
xlim([0 zF])
grid minor
ylabel('Mm, kg/kmol')
xlabel('z, m')

nexttile
plot(z,mu1,'r')
xlim([0 zF])
grid minor
ylabel('\mu , Pa s')
xlabel('z, m')

%% Cinética de reacción
for i=1:size(Fi,1)
r_graf(i,:) = ri(T(i),P(i),Fi(i,:)*1000) ;
R_graf(i,:) = Ri(T(i),P(i),Fi(i,:)*1000) ;
end

figure('Color','White')
t =tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(z,r_graf)
%ylim([0 0.025])
xlim([0 zF])
grid minor
ylabel('r_i , kmol/m^3s')
xlabel('z, m')
legend('Rx. 1)', 'Rx. 2)','Rx. 3)','Rx. 4)','Rx. 5)','Rx. 6)','Rx. 7)')

nexttile
plot(z,R_graf)
xlim([0 zF])
%ylim([-0.001 0.015])
grid minor
ylabel('R_i, kmol/m^3s')
xlabel('z, m')
legend('Metano', 'Acetileno', 'Etileno', 'Etano', 'Propileno', 'Propano', 'Butadieno', 'Hidrogeno', 'Agua')

%% Analisis Perfil Temperatura

figure('Color','White')
t =tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile([2;1])
plot(z,q(z,T)*pi*dt,'r')
%ylim([0 0.025])
xlim([0 zF])
grid minor
ylabel('-q*\piD , kJ/ms')
xlabel('z, m')

nexttile
hold on
for i=1:size(Fi,1)
    FiCp(i) = Fi(i,:)*Cpi' ;
    HR(i) = r_graf(i,:)*Hr(T(i)) ;
end

plot(z,-HR*pi*dt^2/4,'b')
xlim([0 zF])
grid minor
ylabel('-H_R*\piD^2/4, kJ/ms')
xlabel('z, m')

nexttile
plot(z,FiCp,'k')
%ylim([0 0.025])
xlim([0 zF])
grid minor
ylabel('FiCp , kJ/Ks')
xlabel('z, m')

%% RESULTADOS GENERALES
clc

Ntube
F_tubo = F_anual
F_total
zF
Tf
Tg1
Selectividad = Select(end)
Rendimiento = Rend(end)


function cp = Cp(T)

% [  CH4    C2H2     C2H4   C2H6    C3H6     C3H8      C4H8      H2     H2O]
% [Metano Acetileno Etileno Etano Propileno Propano Butadieno Hidrogeno Agua]
% [   1       2        3      4       5        6        7        8       9 ]


C1 = [0.33298, 0.3199, 0.3338, 0.40326, 0.43852, 0.5192, 0.5095, 0.27617, 0.33363]*1e5 ;
C2 = [0.79933, 0.5424, 0.9479, 1.3422, 1.506, 1.9245, 1.705, 0.0956, 0.2679 ]*1e5 ;
C3 = [2.0869, 1.594, 1.596, 1.6555, 1.3988, 1.6265, 1.5324,  2.466, 2.6105 ]*1e3 ;
C4 = [0.41602, 0.4325, 0.551, 0.73223, 0.75754, 1.168, 1.337, 0.0376, 0.08896 ]*1e5 ;
C5 = [991.96, 607.1, 740.8, 752.87, 616.46, 723.6, 685.6, 567.6, 1169] ;

cp = C1 + C2.*((C3./T)./sinh(C3./T)).^2 + C4.*((C5./T)./cosh(C5./T)).^2 ;
cp = cp/1000 ; %kJ/kmolK

end


function mu = Visc(T,Fi)

% [  CH4    C2H2     C2H4   C2H6    C3H6     C3H8     C4H8       H2     H2O]
% [Metano Acetileno Etileno Etano Propileno Propano Butadieno Hidrogeno Agua]
% [   1       2        3      4       5        6       7         8       9 ]

C1 = [5.2546e-7 1.2025e-6 2.0789e-6 2.5906e-7 7.3919e-7 4.9054e-8 2.696e-7 1.797e-7 1.7096e-8] ;
C2 = [0.59006 0.4952 0.4163 0.67988 0.5423 0.90125 0.6715 0.685 1.1146] ;
C3 = [105.67 291.4 352.7 98.902 263.73 0 134.7 -0.59 0] ;
C4 = [0 0 0 0 0 0 0 140 0]  ;

mu = 0 ;

for i=1:9
  mu = mu + Fi(i)/sum(Fi)*(C1(i).*T.^(C2(i))./(1+C3(i)./T + C4(i)./T.^2)) ; %Pa*s
end

%mu = C1.*T.^(C2)./(1+C3./T + C4./T.^2) ;

end


function f = dFidz(i,z,T,P,Fi)
%Balances de materia

global dt

    dFi_dz = Ri(T,P,Fi)*pi*dt^2/4 ; %Vector 1x9
    f = dFi_dz(i) ; %kmol/sm
end


function R = Ri(T,P,Fi)
global Stoi_Coef
    R =  ri(T,P,Fi)*Stoi_Coef ; %vector 1x9 kJ/m3s
end


function r = ri(T,P,Fi)
global R
%Rates ind
% [  CH4    C2H2     C2H4   C2H6    C3H6     C3H8      C4H8      H2     H2O]
% [Metano Acetileno Etileno Etano Propileno Propano Butadieno Hidrogeno Agua]
% [   1       2        3      4       5        6        7        8       9 ]
k1 = @(T) 4.65e13*exp(-273020./(R*T)) ; %s^-1
k2 = @(T) 8.75e8*exp(-136870./(R*T)) ;
k3 = @(T) 3.85e11*exp(-273190./(R*T)) ;
k4 = @(T) 9.81e8*exp(-154580./(R*T)) ;
k5 = @(T) 5.87e4*exp(-29480./(R*T)) ;
k6 = @(T) 1.03e12*exp(-172750./(R*T)) ;
k7 = @(T) 7.08e13*exp(-253010./(R*T)) ;

% PV = NRT -> V = NRT/P

r1 = @(T,P,Fi) k1(T).*Fi(4)./(sum(Fi).*T.*R./P) ;
r2 = @(T,P,Fi) k2(T).*Fi(3).*Fi(8)./(sum(Fi).*T.*R./P).^2 ;
r3 = @(T,P,Fi) k3(T).*Fi(4)./(sum(Fi).*T.*R./P) ;
r4 = @(T,P,Fi) k4(T).*Fi(5)./(sum(Fi).*T.*R./P) ;
r5 = @(T,P,Fi) k5(T).*Fi(2).*Fi(1)./(sum(Fi).*T.*R./P).^2 ;
r6 = @(T,P,Fi) k6(T).*Fi(2).*Fi(3)./(sum(Fi).*T.*R./P).^2 ;
r7 = @(T,P,Fi) k7(T).*Fi(3).*Fi(4)./(sum(Fi).*T.*R./P).^2 ;

r = [r1(T,P,Fi), r2(T,P,Fi), r3(T,P,Fi),...
    r4(T,P,Fi), r5(T,P,Fi), r6(T,P,Fi), r7(T,P,Fi)] ; %Vector 1x7 kmol/sm3

end

