
% SIN REFLUJO
clear
clc

%Datos del equilibrio

%Fase de hidrocarburo

HC_M = [0 9.2 18.6 22 33.8 40.9 46 59 66.8 71.4 ...
    73.3 83.3 88.1]'/100 ; %Metilciclohexano
HC_H = [92.9 83.1 73.4 69.8 57.6 50.4 45 30.7 22.8 ...
    17.8 15.7 5.4 0]'/100 ; %Heptano
HC_A = 1-HC_M-HC_H ; %Anilina

%Fase de hidrocarburo

S_M = [0 0.8 2.7 3 4.6 6 7.4 9.2 11.3 12.7 13.1 ...
    15.6 16.9]'/100; %Metilciclohexano
S_H = [6.2 6 5.3 5.1 4.5 4 3.6 2.8 2.1 1.6 1.4 ...
    0.6 0]'/100; %Heptano
S_A = 1-S_M-S_H ; %Anilina


% Gráfica en coordenadas cartesianas

figure('Color','white')
tiledlayout(1,1)
nexttile

R = [1 0 ; 0 1] ;
%R = [1 cos(pi/3) ; 0 sin(pi/3)] ;

hold on
for i = 0:0.1:1
   Diag = R*[0 i;i 0] ;
   Vert = R*[i i;0 1-i] ;
   Hor = R*[0 1-i;i i];
   
   if i==0 || i==1
   C = [0 0 0] ;
   else
   C = [0.60 0.60 0.60] ;
   end
   
   plot(Diag(1,:),Diag(2,:),'Color',C)
   plot(Vert(1,:),Vert(2,:),'Color',C)
   plot(Hor(1,:),Hor(2,:),'Color',C)
   
   Tx_R = R*[1-i+0.01;i+0.01] ;
   Tx_L = R*[-0.02;1-i+0.01] ;
   Tx_D = R*[i+0.01;-0.03] ;
   
   text(Tx_R(1),Tx_R(2),sprintf('%.0f',i*100))
   text(Tx_L(1),Tx_L(2),sprintf('%.0f',i*100),'HorizontalAlignment','right')
   text(Tx_D(1),Tx_D(2),sprintf('%.0f',i*100))
end

nom = R*[-0.15 -0.05 1.05; 1.05 -0.07 -0.07] ;

text(nom(1,1),nom(2,1),'Metilciclohexano')
text(nom(1,2),nom(2,2),'Heptano')
text(nom(1,3),nom(2,3),'Anilina')
Equi_HC = R*[HC_A';HC_M'] ;
Equi_S = R*[S_A';S_M'] ;

for j = 1:13
    plot([Equi_HC(1,j) Equi_S(1,j)],[Equi_HC(2,j) Equi_S(2,j)],'Color',[0.3 0.3 0.3])
end

plot(Equi_HC(1,:),Equi_HC(2,:),'r',Equi_S(1,:),Equi_S(2,:),'b')
axis('off')

% Composiciones de las corrientes

%Alimentación
wF_H = 0.5 ;
wF_M = 0.5 ;
wF_A = 1 - wF_H - wF_M ;

%Solvente
wS_H = 0.00 ;
wS_M = 0.00 ;
wS_A = 1 - wS_H - wS_M ;

%Refinado
wRNp_H = 0.92 ;
wRNp_M = 0.06 ;
wRNp_A = 1 - wRNp_H - wRNp_M ;

% Transformación de las coordenadas

X = R*[wF_A, wS_A, wRNp_A; wF_M, wS_M, wRNp_M] ;

%Alimentación
wF_A = X(1,1) ;
wF_M = X(2,1) ;
wF_H = 1 - wF_A - wF_M ;

%Solvente
wS_A = X(1,2) ;
wS_M = X(2,2) ;
wS_H = 1 - wS_A - wS_M ;

%Refinado
wRNp_A = X(1,3) ;
wRNp_M = X(2,3) ;
wRNp_H = 1 - wRNp_H - wRNp_M ;

% CÁLCULOS

%Línea entre F y S

LinFS = @(x) interp1([wF_A, wS_A],[wF_M, wS_M],x,'Linear','extrap') ;

%Línea entre R_Np y S

LinRNpS = @(x) interp1([wRNp_A, wS_A],[wRNp_M, wS_M],x,'Linear','extrap') ;

%FUNCIONES DE INTERPOLACIÓN

%Equilibrio M_HC - M_S
S_MEq = @(x) interp1(Equi_HC(2,:),Equi_S(2,:),x,'Linear','extrap') ;
%Equilibrio M_S - A_S
S_A_Eq = @(x) interp1(Equi_S(2,:),Equi_S(1,:),x,'Linear','extrap') ;
%Equilibrio A_HC - M_HC
HC_M_Eq = @(x) interp1(Equi_HC(1,:),Equi_HC(2,:),x,'Linear','extrap') ;
%Equilibrio A_S - M_S
S_M_Eq = @(x) interp1(Equi_S(1,:),Equi_S(2,:),x,'Linear','extrap') ;

% Se crea el vector en el equilbrio HC para extenderlas hasta cruzar con
% la linea R_Np - S
HC_Aeq = linspace(wRNp_M+0.015 ,HC_A(end),20) ;

for i=1:size(HC_Aeq,2)
    LinEq = @(x) interp1([HC_Aeq(i), S_A_Eq(S_MEq(HC_M_Eq(HC_Aeq(i))))],...
        [HC_M_Eq(HC_Aeq(i)), S_MEq(HC_M_Eq(HC_Aeq(i)))],x,'Linear','extrap') ;
    
    obj = @(x) LinEq(x)-LinRNpS(x);
    xinter(i) = fsolve(@(x) LinEq(x)-LinRNpS(x),1) ;
    plot([HC_Aeq(i); xinter(i)], [LinEq(HC_Aeq(i)); LinEq(xinter(i))],'Color',[0.7 0.7 0.7])
end

plot([wRNp_A;wS_A; max(xinter)], [LinRNpS(wRNp_A);LinRNpS(wS_A); LinRNpS(max(xinter))],'g')
axis([-0.1 1.1 -0.1 1.1])

% Se encuentran las coordenadas del delta R min
xDeltaR = min(xinter);
yDeltaR = LinRNpS(xDeltaR);

figure('Color','White')
plot(1:size(xinter,2),xinter)

%% Linea que une el F con el E1_max
LinFEmax = @(x) interp1([wF_A,xDeltaR ],[wF_M, yDeltaR],x,'Linear','extrap') ;
plot([wF_A,xDeltaR ],[wF_M, yDeltaR],'m')

%Coordenadas E1_max
XEmax = fsolve(@(x) LinFEmax(x)-S_M_Eq(x),1);
YEmax = S_M_Eq(XEmax);

%Linea entre R_Np - E1_max
LinLEmax = @(x) interp1([wRNp_A, XEmax ],[wRNp_M, YEmax],x,'Linear','extrap') ;

%Coordenadas punto de mezcla máxima
Mmax_A = fsolve(@(x) LinLEmax(x)-LinFS(x),1);
Mmax_M = LinLEmax(Mmax_A) ;


plot([wRNp_A, XEmax],[wRNp_M, YEmax],'b',Mmax_A,Mmax_M,'o', ...
    [wF_A, wS_A],[wF_M, wS_M],'b')
xlim([0 xDeltaR+0.2])

F=100;
Smin = F*(wF_A-Mmax_A/(Mmax_A-wS_A));

%Intercepciones de las coordenadas de operación desde el punto delta R

for i=1:size(HC_Aeq,2)
    LinEqDelta = @(x) interp1([HC_Aeq(i), xDeltaR],...
        [HC_M_Eq(HC_Aeq(i)), yDeltaR],x,'Linear','extrap') ;
    
 %Anilina del solvente en el punto de operación 
    A_Op(i) = fsolve(@(x) LinEqDelta(x)-S_M_Eq(x),1);
    
end

%Cálculo de S_M en operación con los A calculados
S_Mop = S_M_Eq(A_Op);

%Cálculo de HC_M en operación con los A supuestos para el cálculo
HC_Mop = HC_M_Eq(HC_Aeq);

%Creación de funciones de datos de operación
S_M_Op = @(x) interp1(HC_Mop,S_Mop,x,'Linear','extrap') ;
HC_M_Op = @(x) interp1(S_Mop,HC_Mop,x,'Linear','extrap') ;
HC_MEq = @(x) interp1(S_M,HC_M,x,'Linear','extrap')  ;

%% Etapas

 wE1_M = S_M_Eq(XEmax);

 Xetapas(1) = wF_M ;
 Yetapas(1) = wE1_M ;

 k = 1 ;

 while Xetapas(k) > wRNp_M || k < 100

    Xetapas(k+1) = HC_MEq(Yetapas(k)) ;
    Yetapas(k+1) = S_M_Op(Xetapas(k+1)) ;
    k = k + 1 ;
end

Xetapas = sort([Xetapas,Xetapas]) ;
Yetapas = sort([Yetapas,Yetapas]) ;

nexttile
hold on
plot(wRNp_M:0.01:wF_M,S_M_Op(wRNp_M:0.01:wF_M),'m')
plot([0 1],[0 1],'k',HC_M,S_M,'g')
grid minor
plot(Xetapas(1:end-1),Yetapas(2:end),'k')

