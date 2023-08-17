
%Condiciones TP
P=101.325 ; %KPa
T = 25+273 ;

%Datos
yA1=0.100 ;
xA2=0.005 ;
abs=0.90 ;
m=10.7 ;

%Relaciones molares
XA2=xA2/(1-xA2) ;
YA1=yA1/(1-yA1) ;
YA2=YA1*(1-abs) ;

%Datos calculados
yA2=YA2/(1+YA2) ;
XA1_eq=yA1/(m-yA1) ;
Ls_Gs_min=(YA2-YA1)/(XA2-XA1_eq) ;

%Pendientes de operación.
Ls_Gs1=Ls_Gs_min*1.2 ;
Ls_Gs2=Ls_Gs_min*1.4 ;
Ls_Gs3=Ls_Gs_min*1.6 ;
Ls_Gs4=Ls_Gs_min*1.8 ;

%Composiciones
XA1_1=-(YA2-YA1)/Ls_Gs1+XA2 ;
xA1_1=XA1_1/(1-XA1_1) ;

XA1_2=-(YA2-YA1)/Ls_Gs2+XA2 ;
xA1_2=XA1_2/(1-XA1_2) ;

XA1_3=-(YA2-YA1)/Ls_Gs3+XA2 ;
xA1_3=XA1_3/(1-XA1_3) ;

XA1_4=-(YA2-YA1)/Ls_Gs4+XA2 ;
xA1_4=XA1_4/(1-XA1_4) ;


%Para graficar
X_eq=linspace(0,XA1_eq+XA1_eq/10,41) ;
X_min=linspace(XA2,XA1_eq,41) ;
X_ope1=linspace(XA2,XA1_1,41) ;
X_ope2=linspace(XA2,XA1_2,41) ;
X_ope3=linspace(XA2,XA1_3,41) ;
X_ope4=linspace(XA2,XA1_4,41) ;

x_eq=X_eq./(1+X_eq) ;
Y_eq=m*x_eq./(1-m*x_eq) ;
YA_min=Ls_Gs_min*X_min+YA2 ;

%Curvas de operación

YA_ope1 = Ls_Gs1*X_ope1+YA2 ;
YA_ope2 = Ls_Gs2*X_ope2+YA2 ;
YA_ope3 = Ls_Gs2*X_ope3+YA2 ;
YA_ope4 = Ls_Gs4*X_ope4+YA2 ;


tiledlayout(2,2)

nexttile
plot(X_eq,Y_eq,'k',X_ope1,YA_ope1,'b',X_min,YA_min,'r',...
    XA1_1,YA1,'ob',XA1_eq,YA1,'or',XA2,YA2,'ok')
xlabel('X_{A}') ; ylabel('Y_{A}') ;
legend('Curva Equilibrio','Linea Operación','Min. Operación',...
    'X_{A1},Y_{A1}','X_{Ai},Y_{A1}','X_{A2},Y_{A2}')
grid on
axis tight

nexttile

plot(X_eq,Y_eq,'k',X_ope2,YA_ope2,'b',X_min,YA_min,'r',...
    XA1_2,YA1,'ob',XA1_eq,YA1,'or',XA2,YA2,'ok')
xlabel('X_{A}') ; ylabel('Y_{A}') ;
legend('Curva Equilibrio','Linea Operación','Min. Operación',...
    'X_{A1},Y_{A1}','X_{Ai},Y_{A1}','X_{A2},Y_{A2}')
grid on
axis tight

nexttile

plot(X_eq,Y_eq,'k',X_ope3,YA_ope3,'b',X_min,YA_min,'r',...
    XA1_3,YA1,'ob',XA1_eq,YA1,'or',XA2,YA2,'ok')
xlabel('X_{A}') ; ylabel('Y_{A}') ;
legend('Curva Equilibrio','Linea Operación','Min. Operación',...
    'X_{A1},Y_{A1}','X_{Ai},Y_{A1}','X_{A2},Y_{A2}')
grid on
axis tight

nexttile

plot(X_eq,Y_eq,'k',X_ope4,YA_ope4,'b',X_min,YA_min,'r',...
    XA1_4,YA1,'ob',XA1_eq,YA1,'or',XA2,YA2,'ok')
xlabel('X_{A}') ; ylabel('Y_{A}') ;
legend('Curva Equilibrio','Linea Operación','Min. Operación',...
    'X_{A1},Y_{A1}','X_{Ai},Y_{A1}','X_{A2},Y_{A2}')
grid on
axis tight



%%
FL=1.41e-2 ;
FG=6.84e-4 ;
N=200 ; %Número de nodos integral (Particiones N-1)

yAg=linspace(yA2,yA1,N) ;
factor=1/Ls_Gs*(yAg./(1-yAg)-yA2/(1-yA2)) ;
xAl=factor./(1+factor) ;
yAi=zeros(N,1) ;
F=zeros(N,1) ;

syms x
for i=1:N
     f=((1-x)./(1-yAg(i)))-((1-xAl(i))./(1-x/m)).^(FL/FG) ;
     f=inline(f) ;
     df=inline(diff(f(x))) ;
     yAi0=0.1 ;
    yAi(i)=NR(f,df,yAi0) ;
    F(i)=1/(1-yAg(i))/log((1-yAi(i))/(1-yAg(i))) ;
end
 
result=[yAg' xAl' yAi F] ;
array2table(result,'VariableNames',{'yAg','xAL','yAi','1/(1-y)ln((1-yi)/(1-y))'})
h=(yA1-yA2)/(N-1) ;

   suma=0 ;
   for i=2:N
       suma=suma+F(i) ;
   end
   integral1=h/2*(F(1)+2*suma+F(end)) ;

    Par=0; Impar=0;
    for i=2:N
        if ((-1)^i)==1
           Impar=Impar+F(i) ;
        else
           Par=Par+F(i) ;
        end
     integral2=h/3*(F(1)+4*Impar+2*Par+F(end));
    end
    
 disp('Regla de compuesta simpson')  
sprintf('El area Simpson bajo la curva es: %f ',integral2) 
 disp('Regla de trapecio')
sprintf('El area trapecio bajo la curva es: %f ',integral1)  

tiledlayout(1,2)

nexttile
plot(X_eq,Y_eq,'k',X_ope,YA_ope,'b',X_min,YA_min,'r',...
    XA1,YA1,'ob',XA1_eq,YA1,'or',XA2,YA2,'ok')
xlabel('X_{A}') ; ylabel('Y_{A}') ;
legend('Curva Equilibrio','Linea Operación','Min. Operación',...
    'X_{A1},Y_{A1}','X_{Ai},Y_{A1}','X_{A2},Y_{A2}')
grid on
axis tight

nexttile
plot(yAg',F,'k')
xlabel('y_{A}') ; ylabel({'1';'----------------';'(1-y)ln((1-y_i)/(1-y))'}) ;
hold on
y=[yA2 yAg yA1] ;
A=[0 F' 0] ;
fill(y,A,'b') ;
legend('Area')
grid on
axis tight

