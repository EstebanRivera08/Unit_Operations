clc
clear all

%Condiciones TP

P = [50 100 200 300]; %kPa
%T = [25 35 45 55] ; %°C
T = 25 ; %°C
T = T+273.15 ;
Q = 0.3; %m3/s
R = 8.3145; %kJ/kmol
G = P*Q/R./T*1000; %mol/s

%Nota: 2 hace referencia a la parte de arriba y 1 a la parte de abajo

%Datos
yA1 = 0.100; %Entrada gas
xA2 = 0.005; %Entrada líquido
abs = 0.90; %Fracción absorbida del soluto
m = 4e-5*exp(0.0329*T); %Pendiente de la curva de equilibrio para el pentano

%Relaciones molares

XA2 = xA2/(1-xA2);
YA1 = yA1/(1-yA1);
YA2 = YA1*(1-abs); %Salida del gas

%DATOS CALCULADOS

yA2 = YA2/(1+YA2);

%Composición en el líquido a la salida (flujo mínimo de líquido)
XA1_eq = yA1./(m-yA1); 
xA1_eq = XA1_eq./(1+XA1_eq);

%Pendiente de la curva de operación con flujo mínimo
Ls_Gs_min = (YA2-YA1)./(XA2-XA1_eq);
Gs = G*(1-yA1); %Flujo de gas portador


%Pendientes de operación con cada relación de L respecto a L_min
a= 1.2 ;
Ls_Gs = Ls_Gs_min*a;
Ls = Gs.*Ls_Gs ; %Flujo mínimo de líquido
%Vectores de la curva de equilibrio y de la curva de operación con flujo mínimo
%(Mayúscula-Relacion Molar; Minúscula - Fracción molar)

for i=1:1
X_eq(:,i) = linspace(0,XA1_eq(i)+XA1_eq(i)/10,41); %Valores de x en la curva de equilibrio
X_min(:,i) = linspace(XA2,XA1_eq(i),41); %Valores de X para el flujo mínimo
end

x_eq = X_eq./(1+X_eq) ;
x_min = X_min./(1+X_min);

Y_eq = m.*x_eq./(1-m.*x_eq); %Valores de Y en la curva de equilibrio
YA_min = Ls_Gs_min.*(X_min-XA2)+YA2; %Valores de Y para el flujo mínimo
y_eq = Y_eq./(1+Y_eq);
yA_min = YA_min./(1+YA_min);

%Composición a la salida del líquido (1) con cada flujo 
XA1i = -((YA2-YA1)./Ls_Gs)+XA2; 

%Línea de operación para cada relación de L respecto a L_min
for i=1:1 %filas son temperaturas
    %Valores de X para cada curva de operación
    XA_opei(:,i) = linspace(XA2,XA1i(i),41);
end

%Valores de Y para cada curva de operación
YA_opei = Ls_Gs.*(XA_opei-XA2)+YA2;
%Curvas de operación en fracciones molares
xA1i = XA1i./(1+XA1i);
xA_opei = XA_opei./(1+XA_opei);
yA_opei = YA_opei./(1+YA_opei);


%NÚMERO DE ETAPAS :)
%Matrices con las composiciones de salida de cada corriente luego de cada etapa

X_etapas = zeros(20,4);
X_etapas(1,1:4) = XA2; %Valor de X en la entrada, antes del plato 1

Y_etapas = X_etapas; 
Y_etapas(1,1:4) = YA2; %Valor de Y en la salida, luego del plato 1

y_etapas = X_etapas;
y_etapas(1,1:4) = Y_etapas(1)/(1+Y_etapas(1));

%Cálculo del número de etapas de equilibrio

for j=1:1
    %Composición que se inicia a determinar, ya que XA2 y YA2 ya están determinadas
    i=2; 
    
    %Se calcula una siguiente etapa mientras la concentración de inicio sea
    %menor que XA1
    
    while X_etapas(i-1,j) < XA1i(j) && i<50
                
        %Determinación de X en equilibrio (sale de la etapa)
        X_etapas(i,j) = y_etapas(i-1,j)/(m(j)-y_etapas(i-1,j)); 
        
        %Determinación de Y sobre la curva de operación (entra a la etapa)
        Y_etapas(i,j) = Ls_Gs(j)*(X_etapas(i,j)-XA2)+YA2;
        y_etapas(i,j) = Y_etapas(i,j)/(1+Y_etapas(i,j));
        
        i=i+1 ;
    end
    
    %número de etapas para cada flujo L
    Etapa(j)=find(X_etapas(:,j),1,'last') - 1;
end
x_etapas = X_etapas./(1+X_etapas);

%Matrices con las composiciones de cada etapa en la entrada y la salida duplicadas (para graficar)
X_etapa = [X_etapas(1:Etapa(end)+1,:); X_etapas(1:Etapa(end)+1,:)];
Y_etapa = [Y_etapas(1:Etapa(end)+1,:); Y_etapas(1:Etapa(end)+1,:)];
%Se ordenan de los valores menores a los mayores
X_etapa = sort(X_etapa,1);
Y_etapa = sort(Y_etapa,1);


%GRAFICACIÓN
figure
t=tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';


for j=1:4
    nexttile
    %Posición en X_Etapa del primer elemento de concentración
    inicio = find(X_etapa(:,1),1);
    %Posición en X_Etapa del último elemento de concentración
    final = find(X_etapa(:,1),1,'last');
    
    %Curvas y vectores a graficar    
    plot(X_eq(:,1),Y_eq(:,1),'k',XA_opei(:,1),YA_opei(:,1),'b',X_min(:,1),YA_min(:,1),'r',...
    XA1i(1),YA1,'ob',XA1_eq(1),YA1,'or',XA2,YA2,'ok',...
    X_etapa(inicio+1:final,1),Y_etapa(inicio:final-1,1),'k')
    
    xlabel('X_{A}'); 
    ylabel('Y_{A}');
    
    grid on
    axis tight
    
   if j==1 
    legend({'Curva Equilibrio','Linea Operación','Min. Operación',...
    'X_{A_{n+1}},Y_{A_{0}}','X_{Ai},Y_{A1}','X_{A_{0}},Y_{A_{n+1}}'},'Location','southeast')
   end
   txt = {sprintf('No. de etapas: %d',Etapa(1)),...
       sprintf('Flujo líquido: %0.2f mol/s',Ls(j)),sprintf('Flujo gas: %0.2f mol/s',Gs(j))} ;
   text(0.005,0.1,txt) 
   title(sprintf('P_{%d} = %0.0f kPa',j,P(j)))
end
    
    %Gráfica del perfil de concentraciones
figure
t=tiledlayout(2,2) ;
t.TileSpacing = 'compact';
t.Padding = 'compact';

%Se grafica la etapa (Eta) y la concentración del gas que entra y del 
%líquido que sale de esta

%Nota: Eta 0 hace referencia a las condiciones del gas que sale y del líq
%que entra a la columna

for j=1:4
    nexttile

    Eta = find(X_etapas(:,1))-1; %Etapas a graficar para cada caso
    
    plot(Eta,X_etapas(Eta+1,1),Eta,Y_etapas(Eta+1,1),Eta,ones(size(Eta))*XA1i(1),'--',Eta,ones(size(Eta))*YA1,'--') 

    xlabel('Etapa'); 
    ylabel('Relación molar');
    
    if j==1
    legend({'X_{A}','Y_{A}','X_{A0}','Y_{A0}'}, 'Location','southeast')
    end
    title(sprintf('P_{%d} = %0.0f kPa',j,P(j)))
    grid on
    axis tight

end   
