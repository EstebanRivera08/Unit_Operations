clear all 
clc

%Componentes
%1 - Etano
%2 - Propano
%3 - n-Butano
%4 - n-Pentano
%5 - n-Hexano

F = zeros(16, 1);
W = zeros(16, 1);
V = zeros(17, 1);
L = zeros(16, 1);
Q = zeros(16, 1);
U = zeros(16, 1);
global T P x
T = zeros(16, 1);
P = zeros(16, 1);
x = zeros(16, 5);
y = zeros(17, 5);
h = zeros(16, 1);
H = zeros(16, 1);
zF = zeros(16, 5); 
TF = zeros(16, 1); 
PF = zeros(16, 1); 
hF = zeros(16, 1);

%Valores definidos --------------------------------------------------------

%Flujos de entrada -- lbmol/h
F(6) = 41; F(9) = 59;

%Flujos de entrada -- Composiciones
zF(6, 1) = 2.5/41; zF(6, 2) = 14/41; zF(6, 3) = 19/41; zF(6, 4) = 5/41; zF(6, 5) = 0.5/41;
zF(9, 1) = 0.5/59; zF(9, 2) = 6/59; zF(9, 3) = 18/59; zF(9, 4) = 30/59; zF(9, 5) = 4.5/59;

%Flujos de entrada -- Temperatura -- °F
TF(6) = 170; TF(9) = 230;

%Flujos de entrada -- Presión -- psia
PF(6) = 300; PF(9) = 275;

%Presiones -- psia
P(1) = 238; P(16) = 243;
P(2) = 240; 
for i=3:1:15
    P(i) = 240 + (i-2)*0.2;
end

%Productos liq intermedios -- lbmol/h
U(1) = 5; U(3) = 3;

%Productos vap intermedios -- lbmol/h
W(13) = 37;

%Flujos de calor -- btu/h
Q(3) = 200000;

%Vapor salida de 1 lbmol/h
V(1) = 15;

%Liquido salida de 1 a 2
L(1) = 150;


%SUPOSICIONES INICIALES ---------------------------------------------------

%Sup - Salida superior
x(1,1)=3/42; x(1,2)=17/42; x(1,3)=20/42; x(1,4)=2/42; x(1,5)=0/42;

%Sup - Salida inferior
x(16,1)=0/58; x(16,2)=3/58; x(16,3)=17/58; x(16,4)=33/58; x(16,5)=5/58;

%Determinación de temperatura
Kk = @(T,n) K(1,P(n),T)*x(n,1)+K(2,P(n),T)*x(n,2)+K(3,P(n),T)*x(n,3)+...
    K(4,P(n),T)*x(n,4)+K(5,P(n),T)*x(n,5)-1;

Tmin = fsolve(@(T) Kk(T,1), 50);
Tmax = fsolve(@(T) Kk(T,16), 250);

%Temperaturas -- °F
for i=1:1:16
    T(i) = Tmin + (i - 1)*((Tmax-Tmin)/15);
end

%Flujos de vapor a etapa superior - lbmol/h
V(2) = L(1) + U(1) + V(1);
for i=3:1:6
    V(i) = V(2);
end
Kf = @(T,n) K(1,PF(n),T)*zF(n,1)+K(2,PF(n),T)*zF(n,2)+K(3,PF(n),T)*zF(n,3)+...
    K(4,PF(n),T)*zF(n,4)+K(5,PF(n),T)*zF(n,5)-1;

TF6eb = fsolve(@(T) Kf(T,6), TF(6));
TF9eb = fsolve(@(T) Kf(T,9), TF(9));

if TF6eb<TF(6)
    V(7) = V(6) - F(6);
    hF(6) = 0;
    for i=1:1:5
        hF(6) = hF(6)+ zF(6,i)*(lambda(i,0) + integral(@(Tn) Cpvap(i,Tn),0,TF(6)));
    end
else
    V(7) = V(6);
    hF(6)=0;
    for i=1:1:5
        hF(6) = hF(6) + zF(6,i)*Cpliq(i, TF(6)/2)*TF(6);
    end
end
V(8) = V(7); V(9) = V(7);
if TF9eb<TF(9)
    V(10) = V(9) - F(9);
    hF(9) = 0;
    for i=1:1:5
        hF(9) = hF(9)+ zF(9,i)*(lambda(i,0) + integral(@(Tn) Cpvap(i,Tn),0,TF(9)));
    end
else
    V(10) = V(9);
    hF(9)=0;
    for i=1:1:5
        hF(9) = hF(9) + zF(9,i)*Cpliq(i, TF(9)/2)*TF(9);
    end
end
V(11) = V(10); V(12) = V(10);
V(13) = V(12) + W(13);
V(14) = V(13); V(15) = V(13); V(16) = V(13);


%INICIO DE ITERACIÓN-------------------------------------------------------
tao=10;
iter = 0;
Temperaturas = T ;
ComposicionesLiq(:,:,1) = x ;

while tao > 0.00001
    
    iter = iter + 1 ;
    
    Tanterior = T ; 
    
    %MATRIZ TRIDIAGONAL--------------------------
    
    %Vector de A
    A = zeros(16,1);
    for j=2:1:16
        SumA(j) = 0;
        for m=1:1:j-1
            SumA(j)=SumA(j) + F(m)-U(m)-W(m);
        end
        A(j) = V(j) + SumA(j) - V(1);
    end
    
    %Vector de B
    B = zeros(16,5);
    for i=1:1:5
        for j=1:1:16
            SumB(j)=0;
            for m=1:1:j
                SumB(j)= SumB(j) + F(m)-U(m)-W(m);
            end
            B(j,i) = -( V(j+1) + SumB(j) - V(1) + U(j) + (V(j)+W(j))*K(i,P(j),T(j)) );
        end
    end
    
    %Vector de C
    C = zeros(16,5);
    for i=1:1:5
        for j=1:1:15
            C(j,i) = V(j+1)*K(i,P(j+1),T(j+1));
        end
    end
    
    %Vector de D
    D = zeros(16,5);
    for i=1:1:5
        for j=1:1:16
            D(j,i) = -F(j)*zF(j,i);
        end
    end

    %Matrices tridiagonales
    Tri1 = zeros(16); Tri2 = zeros(16); Tri3 = zeros(16); Tri4 = zeros(16); Tri5 = zeros(16);
    
    for f=1:1:16
        for c=1:1:16
            if f == c
                Tri1(f,c)=B(f,1);
                Tri2(f,c)=B(f,2);
                Tri3(f,c)=B(f,3);
                Tri4(f,c)=B(f,4);
                Tri5(f,c)=B(f,5);
            end
            if f == c+1
                Tri1(f,c)=A(f);
                Tri2(f,c)=A(f);
                Tri3(f,c)=A(f);
                Tri4(f,c)=A(f);
                Tri5(f,c)=A(f);
            end
            if f == c-1
                Tri1(f,c)=C(f,1);
                Tri2(f,c)=C(f,2);
                Tri3(f,c)=C(f,3);
                Tri4(f,c)=C(f,4);
                Tri5(f,c)=C(f,5);
            end
        end
    end
    
    D1 = D(:,1);
    x1 = inv(Tri1)*D1;
    D2 = D(:,2);
    x2 = inv(Tri2)*D2;
    D3 = D(:,3);
    x3 = inv(Tri3)*D3;
    D4 = D(:,4);
    x4 = inv(Tri4)*D4;
    D5 = D(:,5);
    x5 = inv(Tri5)*D5;
    
    x = [x1 x2 x3 x4 x5];

    %Normalizar composiciones---------------------
    Sum = sum(x,2);
    for j=1:1:16
        for i=1:1:5
            x(j,i)=x(j,i)./Sum(j);
        end    
    end
    clc
    %Nuevas temperaturas--------------------------
    for j=1:1:16
        T(j) = fsolve(@(T) SolveT(j,T), T(j));
    end
    
    %Nuevos y-------------------------------------
    for j=1:1:16
        for i=1:1:5
            y(j,i) = K(i, P(j), T(j))*x(j,i);
        end
    end
    
    %Calores en el rehervidor y condensador--------------
    
    for j=1:1:16
        h(j) = hliq(j);
        H(j) = hvap(j);
    end
    
    Q(1) = V(2)*H(2) - (L(1)+U(1))*h(1) - V(1)*H(1); %Condensador
    
    sumQ = 0;
    for j=1:16
        sumQ = sumQ + F(j)*hF(j) - U(j)*h(j) - W(j)*H(j);
    end
    sumQ2 = 0;
    for j=1:15
        sumQ2 = sumQ2 + Q(j);
    end
    
    Q(16) = sumQ - sumQ2 - V(1)*H(1) - L(16)*h(16);
    
    %Nuevos V--------------------------------
    
    %Vector de Alpha
    Alpha = zeros(16,1);
    Alpha(1) = -H(1);
    for j=2:1:16
        Alpha(j) = h(j-1)-H(j);
    end
    
    %Vector de beta
    beta = zeros(16,1);
    for j=1:1:15
        beta(j) = H(j+1)-h(j);
    end
    
    %Vector gamma
    gamma = zeros(16,1);
    Sumgamma = zeros(16,1);
    
    Summgama(1) = 0;
    gamma(1) = (Sumgamma(1) - V(1))*(h(1)-0) + F(1)*(h(1)-hF(1))...
            + W(1)*(H(1)-h(1)) + Q(1);
    
    for j=2:1:16
        for m=1:1:j-1
            Sumgamma(j)= Sumgamma(j) + F(m)-U(m)-W(m);
        end
        gamma(j) = (Sumgamma(j) - V(1))*(h(j)-h(j-1)) + F(j)*(h(j)-hF(j))...
            + W(j)*(H(j)-h(j)) + Q(j);
    end

    for j=2:16
        V(j) = (gamma(j-1) - V(j-1)*Alpha(j-1))/beta(j-1);
    end
    
    %Nuevos L----------------------------------
    SumL=zeros(16,1);
    for j=2:16
        SumL(j)=0;
        for m=1:1:j
            SumL(j)= SumL(j) + F(m)-U(m)-W(m);
        end
        L(j) = V(j+1)+SumL(j) - V(1);
    end
    
    %Verificación-----------------------------
    
    tao = 0;
    for j=1:16
        tao = tao + ((T(j)-Tanterior(j))^2);
    end
    Diferencia(iter) = tao ;
    Temperaturas = [Temperaturas T] ;
    ComposicionesLiq(:,:,iter+1) = x ; 
    ComposicionesVap(:,:,iter) = y ; 
end


%% GRÁFICOS ----------------------------------------------------------
%1 - Etano
%2 - Propano
%3 - n-Butano
%4 - n-Pentano
%5 - n-Hexano

%RESUMEN DE GRÁFICAS

figure('Color','White')
Etapas = 1:16 ;
t = tiledlayout(2,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(Etapas',x)
legend('Etano','Propano','n-Butano','n-Pentano','n-Hexano')
xlabel('No. Etapas')
ylabel('Composiciones, x_i')
xlim([1 16])
grid minor

nexttile
plot(Etapas',y(1:16,:))
legend('Etano','Propano','n-Butano','n-Pentano','n-Hexano')
xlabel('No. Etapas')
ylabel('Composiciones, y_i')
xlim([1 16])
grid minor

nexttile
plot(Etapas',T)
xlabel('No. Etapas')
ylabel('Temperatura, °F')
xlim([1 16])
grid minor

nexttile
plot(Etapas',[L V(1:16)],Etapas',[F -W -U])
legend('Liquido', 'Vapor','F','W','U')
xlabel('No. Etapas')
ylabel('Flujo molar, lbmol/h')
xlim([1 16])
grid minor

%%

figure('Color','White')
Etapas = 1:16 ;
t = tiledlayout(2,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile([2 1])
plot(Etapas',[L.*h V(1:16).*H])
legend('L*h_{Liquido}', 'V*H_{Vapor}')
xlabel('No. Etapas')
ylabel('Entalpía, BTU/h')
xlim([1 16])
grid minor

nexttile
plot(Etapas',[h H])
legend('h_{Liquido}', 'H_{Vapor}')
xlabel('No. Etapas')
ylabel('Entalpía, BTU/lbmolh')
xlim([1 16])
grid minor

nexttile
plot(Etapas',Q,'k')
xlabel('No. Etapas')
ylabel('Q, BTU/h')
xlim([1 16])
grid minor

%%
%COMPARACIONES ENTRE ITERACIONES

Iteraciones = 1:iter ;

%Tau vs. iteraciones
figure('Color','White')
plot(Iteraciones,Diferencia,'k')
xlabel('Iteraciones')
ylabel('\tau = suma[ T_{i}-T_{i-1} ]')
grid on

%T vs iteraciones
figure('Color','White')
t = tiledlayout(1,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(Etapas',T,'k')
xlabel('No. Etapas')
ylabel('Temperaturas, °F')
xlim([1 16])
grid minor

nexttile
for j = 1:6
hold on
plot(Etapas',Temperaturas(:,j))
end
xlabel('No. Etapas')
ylabel('Temperaturas, °F')
grid minor
legend('i = 0','i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off

%x vs iteraciones

figure('Color','White')
t = tiledlayout(3,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(Etapas',x)
legend('Etano','Propano','n-Butano','n-Pentano','n-Hexano')
xlabel('No. Etapas')
ylabel('Composiciones, x_i')
xlim([1 16])
grid minor

nexttile
for j = 1:6
hold on
plot(Etapas',ComposicionesLiq(:,1,j))
end
xlabel('No. Etapas')
ylabel('x_{Etano}')
grid minor
legend('i = 0','i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off

nexttile
for j = 1:6
hold on
plot(Etapas',ComposicionesLiq(:,2,j))
end
xlabel('No. Etapas')
ylabel('x_{Propano}')
grid minor
legend('i = 0','i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off

nexttile
for j = 1:6
hold on
plot(Etapas',ComposicionesLiq(:,3,j))
end
xlabel('No. Etapas')
ylabel('x_{n - Butano}')
grid minor
legend('i = 0','i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off

nexttile
for j = 1:6
hold on
plot(Etapas',ComposicionesLiq(:,4,j))
end
xlabel('No. Etapas')
ylabel('x_{n - Pentano}')
grid minor
legend('i = 0','i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off

nexttile
for j = 1:6
hold on
plot(Etapas',ComposicionesLiq(:,5,j))
end
xlabel('No. Etapas')
ylabel('x_{n - Hexano}')
grid minor
legend('i = 0','i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off

% y vs iteraciones

figure('Color','White')
t = tiledlayout(3,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(Etapas',y(1:16,:))
legend('Etano','Propano','n-Butano','n-Pentano','n-Hexano')
xlabel('No. Etapas')
ylabel('Composiciones, y_i')
xlim([1 16])
grid minor

nexttile
for j = 1:5
hold on
plot(Etapas',ComposicionesVap(1:16,1,j))
end
xlabel('No. Etapas')
ylabel('y_{Etano}')
grid minor
legend('i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off

nexttile
for j = 1:5
hold on
plot(Etapas',ComposicionesVap(1:16,2,j))
end
xlabel('No. Etapas')
ylabel('y_{Propano}')
grid minor
legend('i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off

nexttile
for j = 1:5
hold on
plot(Etapas',ComposicionesVap(1:16,3,j))
end
xlabel('No. Etapas')
ylabel('y_{n - Butano}')
grid minor
legend('i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off

nexttile
for j = 1:5
hold on
plot(Etapas',ComposicionesVap(1:16,4,j))
end
xlabel('No. Etapas')
ylabel('y_{n - Pentano}')
grid minor
legend('i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off

nexttile
for j = 1:5
hold on
plot(Etapas',ComposicionesVap(1:16,5,j))
end
xlabel('No. Etapas')
ylabel('y_{n - Hexano}')
grid minor
legend('i = 1','i = 2','i = 3','i = 4','i = 5')
xlim([1 16])
hold off


%%


%Verificación balances materia
Ver  = zeros(16,5);
for i=1:5
    Ver(1,i) = L(1)*x(1,i) + V(1)*y(1,i) + U(1)*x(1,i) + W(1)*y(1,i) - V(2)*y(2,i);
end

for j=2:16
    for i=1:5
        Ver(j,i) = L(j)*x(j,i) + V(j)*y(j,i) + U(j)*x(j,i) + W(j)*y(j,i)...
            - F(j)*zF(j,i) - L(j-1)*x(j-1,i) - V(j+1)*y(j+1,i);
    end
end

%FUNCIONES-----------------------------------------------------------------

%Entalpía líquido saturado
function entalpia = hliq(n)
    global T x
    entalpia=0;
    for i=1:1:5
        entalpia = entalpia + x(n,i)*Cpliq(i, T(n)/2)*T(n);
    end
end

%Entalpía vapor saturado
function entalpia = hvap(n)
    global T x
    entalpia=0;
    for i=1:1:5
        entalpia = entalpia + x(n,i)*(lambda(i,0) + integral(@(Tn) Cpvap(i,Tn),0,T(n)));
    end
end

%Cps líquidos de los componentes
function cepe = Cpliq(i, Tv)
    
    if i == 1 && Tv > 62.33
        Tv = 62.33;
    end
    if i == 2 && Tv > 188.33
        Tv = 188.33;
    end
    if i == 3 && Tv > 260.33
        Tv = 260.33;
    end
    if i == 4 && Tv > 242.33
        Tv = 242.33;
    end
    if i == 5 && Tv > 368.33
        Tv = 368.33;
    end
    
    TR = Tv + 459.67;
    TK = TR/1.8; 
    TcK = [305.32 369.83 425.12 469.7 507.6];
    TcR = 1.8*TcK;
    Tr=TR./TcR;
    t = 1-Tr;
    C1 = [44.009 62.983 191030 159080 172120];
    C2 = [89718 113630 -1675 -270.5 -183.78];
    C3 = [918.77 633.21 12.5 9.9537e-1 8.8734e-1];
    C4 = [-1886 -873.46 -0.03874 0 0];
    C5 = [0 0 4.6121e-5 0 0];
    if i == 1 || i == 2
        cepe1 = C1(i).*C1(i)./t(i) + C2(i) - (2.*C1(i).*C3(i)).*t(i)...
            - C1(i).*C4(i).*(t(i).^2) - ((C3(i).^2).*(t(i).^3))/3 -...
            (t(i).^4).*C3(i).*C4(i)./2 - ((C4(i).^2).*(t(i).^5))/5;
    elseif i == 4 || i == 5 || i == 3
        cepe1 = C1(i) + C2(i).*TK + C3(i).*TK.*TK + C4(i).*TK.*TK.*TK + C5(i).*(TK.^4);
    end
    cepe = cepe1*2.390058e-4; %btu/lbmol*F
end

%Cps vapor de los componentes
function cepe = Cpvap(i, T)
    TR = T + 459.67;
    TK = TR/1.8;
    
    C1 = [0.4033 0.5192 0.7134 0.8805 1.0440].*1e5;
    C2 = [1.3422 1.9245 2.43 3.011 3.523].*1e5;
    C3 = [1.6555 1.6265 1.63 1.6502 1.6946].*1e3;
    C4 = [0.7322 1.168 1.5033 1.8920 2.369].*1e5;
    C5 = [752.87 723.6 730.42 747.6 761.6];
    
    cepe1 = C1(i) + C2(i).*((C3(i)./(TK.*sinh(C3(i)./TK))).^2) +...
        C4(i).*((C5(i)./(TK.*cosh(C5(i)./TK))).^2);
    
    cepe = cepe1*2.390058e-4; %btu/lbmol*F
end

%Calor de vaporización 
function f=lambda(i, T) 
    C1 = [2.1091 2.9209 3.6238 3.9109 4.4544].*1e7;
    C2 = [0.60646 0.78237 0.8337 0.38681 0.39002];
    C3 = [-0.55492 -0.77319 -0.82274 0 0];
    C4 = [0.32799 0.39246 0.39613 0 0];
    TR = T + 459.67;
    Tr = [TR/305.32 TR/369.83 TR/425.12 TR/469.7 TR/507.6]./1.8;
    
    f = C1(i)*((1-Tr(i))^(C2(i)+C3(i)*Tr(i)+C4(i)*Tr(i)*Tr(i)))*4.302106e-4; %btu/lbmol
end

%Cálculo de la temperatura en una etapa
function SolveTemp = SolveT(n, T)
global x;
global P;
    SolveTemp = K(1,P(n),T)*x(n,1)+K(2,P(n),T)*x(n,2)+K(3,P(n),T)*x(n,3)+...
    K(4,P(n),T)*x(n,4)+K(5,P(n),T)*x(n,5)-1;
end

%Coeficientes de distribuición McWilliams
function KMcWilliams = K(i, P, Tf)

at1 = [-687248.25 -970688.5625 0 -1644864 0];
at2 = [0 0 0 0 0];
at3 = [0 0 19.65479 0 0];
at4 = [0 0 -0.02024 0 0.04476];
at5 = [0 0 0 0 -0.0000233488];
at6 = [7.90699 7.71725 -109.11067 8.3288 -15.52781];
ap1 = [-0.886 -0.67984 -0.99838 -1.17078 -1.23197];
ap2 = [49.02654 0 0 0 0];
ap3 = [0 6.90224 0 0 0];
ap4 = [0 0 0 0 0];
ap5 = [0 0 0 0.00523 0.00718];
ap6 = [0 0 0 0 0];

%P-psia
T = Tf + 459.67; %R

KMcWilliams = exp( at1(i).*(1./(T.^2)) + at2(i).*(1./(T)) +at3(i).*log(T) + ...
    at4(i).*(T) + at5(i).*(T.^2) + at6(i) + ap1(i).*log(P) + ap2(i).*(1./(P.^2)) ...
    + ap3(i).*(1./(P)) + ap4(i).*((log(P)).^2) + + ap5(i).*((log(P)).^3) + ap6(i).*P );
end


