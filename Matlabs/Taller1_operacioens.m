clc
clear all

%Condiciones TP
P = 101.325*860/760 ; %kPa
T = 25+273 ; %K
Q = 0.3 ; %m3/s
R = 8.3145 ; %kJ/kmol
G = P*Q/R/T*1000 ; %mol

%Datos
yA1=0.100 ;
xA2=0.005 ;
abs=0.90 ;
m=4e-5*exp(0.0329*T) ;

%Relaciones molares
XA2=xA2/(1-xA2) ;
YA1=yA1/(1-yA1) ;
YA2=YA1*(1-abs) ;

%Datos calculados
yA2=YA2/(1+YA2) ;
XA1_eq=yA1/(m-yA1) ;
xA1_eq=XA1_eq/(1+XA1_eq) ;
Ls_Gs_min=(YA2-YA1)/(XA2-XA1_eq) ;
Gs = G*(1-yA1) ;
Ls_min = Gs*Ls_Gs_min ;

%Pendientes de operación.
Ls_Gs(1,1) = Ls_Gs_min*1.2 ;
Ls_Gs(1,2) = Ls_Gs_min*1.4 ;
Ls_Gs(1,3) = Ls_Gs_min*1.6 ;
Ls_Gs(1,4)=Ls_Gs_min*1.8 ;

%VALORES PARA GRAFICAR

%Relaciones
X_eq=linspace(0,XA1_eq+XA1_eq/10,41) ;
X_min=linspace(XA2,XA1_eq,41) ;

YA_min=Ls_Gs_min*(X_min-XA2)+YA2 ;
%Composiciones
x_eq = X_eq./(1+X_eq) ;
x_min = X_min./(1+X_min) ;
Y_eq = m*x_eq./(1-m*x_eq) ;
y_eq = Y_eq./(1+Y_eq) ;
yA_min = YA_min./(1+YA_min) ;

for i=1:4
%curva de operación en relaciones
XA1i(1,i)=-(YA2-YA1)/Ls_Gs(i)+XA2 ;
XA_opei(:,i)=linspace(XA2,XA1i(1,i),41) ;
YA_opei(:,i) = Ls_Gs(i)*(XA_opei(:,i)-XA2)+YA2 ;
end

%curva de operación en relaciones
xA1i = XA1i./(1+XA1i) ;
xA_opei = XA_opei./(1+XA_opei) ;
yA_opei = YA_opei./(1+YA_opei) ;

%Número de etapas
X_etapas = zeros(20,4) ;
Y_etapas = X_etapas ; y_etapas = X_etapas ;

X_etapas(1,1:4) = XA2 ;
Y_etapas(1,1:4) = YA2 ; 
y_etapas(1,1:4) = Y_etapas(1)/(1+Y_etapas(1)) ;

for j=1:4
    i=2 ;
    while X_etapas(i-1,j) < XA1i(j) && i<20
       y_etapas(i,j) = Y_etapas(i-1,j)/(1+Y_etapas(i-1,j)) ;
       X_etapas(i,j) = y_etapas(i,j)/(m-y_etapas(i,j)) ;
       Y_etapas(i,j) = Ls_Gs(j)*(X_etapas(i,j)-XA2)+YA2 ;
    i=i+1 ;
    end
Etapa(j) = find(X_etapas(:,j), 1, 'last')-1 ;
end
x_etapas = X_etapas./(1+X_etapas) ;
X_etapa = [X_etapas(1:Etapa(1)+1,:) ; X_etapas(1:Etapa(1)+1,:)] ;
Y_etapa = [Y_etapas(1:Etapa(1)+1,:) ; Y_etapas(1:Etapa(1)+1,:)] ;

X_etapa = sort(X_etapa,1) ;
Y_etapa = sort(Y_etapa,1) ;

figure
t=tiledlayout(2,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';


for j=1:4 
  nexttile
inicio = find(X_etapa(:,j),1) ;
final = find(X_etapa(:,j),1,'last') ;

plot(X_eq,Y_eq,'k',XA_opei(:,j),YA_opei(:,j),'b',X_min,YA_min,'r',...
     X_etapa(inicio+1:final,j),Y_etapa(inicio:final-1,j),'k',...
        XA1i(j),YA1,'ob',XA1_eq,YA1,'or',XA2,YA2,'ok')
      
xlabel('X_{A}') ; ylabel('Y_{A}') ;
grid on
axis tight
    if j==1 
    legend('Curva Equilibrio','Linea Operación','Min. Operación',...
    'X_{A1},Y_{A1}','X_{Ai},Y_{A1}','X_{A2},Y_{A2}')
    title('a = 1.2')
    elseif j==2 
    title('a = 1.4')
    elseif j==3
    title('a = 1.6')
    elseif j==4
    title('a = 1.8')
    end
end   

figure
t=tiledlayout(2,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';


for j=1:4 
  nexttile

Eta = find(X_etapas(:,j))-1 ;
plot(Eta,X_etapas(Eta+1,j),'m',Eta,Y_etapas(Eta+1,j),'g',...
    Eta,ones(size(Eta))*XA1i(j),'m--',Eta,ones(size(Eta))*YA1,'g--')

 xlabel('No. Etapa') ; ylabel('X_{A},Y_{A}') ;

grid on
axis tight

    if j==1 
    legend('X_{A}','Y_{A}','X_{A1}','Y_{A1}')
    title('a = 1.2')
    elseif j==2 
    title('a = 1.4')
    elseif j==3
    title('a = 1.6')
    elseif j==4
    title('a = 1.8')
    end

end   

fprintf('El flujo mínimo de aceite es: %f \n',Ls_min)
%xy etampas




