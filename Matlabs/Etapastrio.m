

F = 10 ;
E = 10 ;
R = 10 ;
S = 10 ;

load('TernarySystemS2.mat')
HC_Metilciclohexano = table2array(TernarySystemS2(1:13,1)) ;
HC_Heptano = table2array(TernarySystemS2(1:13,2)) ;
HC_Anilina = table2array(TernarySystemS2(1:13,3)) ;

S_Metilciclohexano = table2array(TernarySystemS2(1:13,5)) ;
S_Heptano = table2array(TernarySystemS2(1:13,6)) ;
S_Anilina = table2array(TernarySystemS2(1:13,7)) ;


figure('Color','white')

hold on
for i = 0:10:100
   Diag = [0 i;i 0] ;
   Vert = [i i;0 100-i] ;
   Hor = [0 100-i;i i];
   
   if i==0 || i==100
   C = [0 0 0] ;
   else
   C = [0.60 0.60 0.60] ;
   end
   
   plot(Diag(1,:),Diag(2,:),'Color',C)
   plot(Vert(1,:),Vert(2,:),'Color',C)
   plot(Hor(1,:),Hor(2,:),'Color',C)
   
   Tx_R = [100-i+1;i+1] ;
   Tx_L = [-2;100-i+1] ;
   Tx_D = [i+1;-3] ;
   
   text(Tx_R(1),Tx_R(2),sprintf('%.0f',i))
   text(Tx_L(1),Tx_L(2),sprintf('%.0f',i),'HorizontalAlignment','right')
   text(Tx_D(1),Tx_D(2),sprintf('%.0f',i))
end

for j = 1:13
    plot([HC_Anilina(j) S_Anilina(j)],[HC_Metilciclohexano(j) S_Metilciclohexano(j)],'Color',[0.3 0.3 0.3])
end

nom = [-10 -10 100; 105 -7 -7] ;

text(nom(1,1),nom(2,1),'Metilciclohexano')
text(nom(1,2),nom(2,2),'Heptano')
text(nom(1,3),nom(2,3),'Anilina')
plot(HC_Anilina,HC_Metilciclohexano,'r',S_Anilina,S_Metilciclohexano,'b')
axis('off')

xF_A = 0 ;
xF_H = 0.5 ;
xF_M = 1 - xF_H - xF_A ;

xS_A = 1 ;
xS_H = 0 ;
xS_M = 1 - xS_H - xS_A ;

xE_M = 0.01 ;
%BM
M = S + F ;

xM_A = (S*xS_A + F*xF_A)/(S + F) ;
xM_M = (S*xS_M + F*xF_M)/(S + F) ;
xM_H = 1 - xM_A - xM_M ;
xRN_M = 0.01 ;

xS_M_eq = @(x) interp1(S_Anilina/100,S_Metilciclohexano/100,x,'linear','extrap') ;
xR_M_eq = @(x) interp1(HC_Anilina/100,HC_Metilciclohexano/100,x,'linear','extrap') ; 
xS_A_eq = @(x) interp1(S_Metilciclohexano/100,S_Anilina/100,x,'linear','extrap') ;
xR_A_eq = @(x) interp1(HC_Metilciclohexano/100,HC_Anilina/100,x,'linear','extrap') ; 

x_A = @(x) interp1(S_Anilina/100,HC_Anilina/100,x,'linear','extrap') ;

xAnil = linspace(min(S_Anilina)/100,max(S_Anilina)/100,20) ;

for i = 1:size(xAnil,2) 
LinEq = @(x) (xS_M_eq(xAnil(i)) - xR_M_eq(x_A(xAnil(i))))./(xAnil(i)-...
    x_A(xAnil(i)))*(x-xAnil(i)) + xS_M_eq(xAnil(i)) ;

LinRS = @(x) (xS_M-xRN_M)/(xS_A-xR_A_eq(xRN_M))*(x - xR_A_eq(xRN_M)) + xRN_M ;

deltaR(i) = fsolve(@(x) LinEq(x)-LinRS(x), 1) ;

plot([xAnil(i) x_A(xAnil(i)) deltaR(i)],[xS_M_eq(xAnil(i))...
    xR_M_eq(x_A(xAnil(i))) LinRS(deltaR(i))],'o-g')

plot([xS_A xM_A xF_A]*100,[xS_M xM_M xF_M]*100,'o-g')

end



hold on
plot([xS_A xM_A xF_A]*100,[xS_M xM_M xF_M]*100,'o-g')
%plot([xF_M_eq xM_A xF_A]*100,[xS_M xM_M xF_M]*100,'o-g')


