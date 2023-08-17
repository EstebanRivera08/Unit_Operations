clc
clear

load('TernarySystemS2.mat')

HC_Metilciclohexano = table2array(TernarySystemS2(1:13,1)) ;
HC_Heptano = table2array(TernarySystemS2(1:13,2)) ;
HC_Anilina = table2array(TernarySystemS2(1:13,3)) ;

S_Metilciclohexano = table2array(TernarySystemS2(1:13,5)) ;
S_Heptano = table2array(TernarySystemS2(1:13,6)) ;
S_Anilina = table2array(TernarySystemS2(1:13,7)) ;

figure('Color','white')
tiledlayout(1,2)

nexttile
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

nexttile
R = [1 cos(pi/3) ; 0 sin(pi/3)] ;

hold on
for i = 0:10:100
   Diag = [0 i;i 0] ;
   Vert = [i i;0 100-i] ;
   Hor = [0 100-i;i i];
   
   Diag = R*Diag ;
   Vert = R*Vert;
   Hor =  R*Hor;
   
   if i==0 || i==100
   C = [0 0 0] ;
   else
   C = [0.50 0.50 0.50] ;
   end
   
   plot(Diag(1,:),Diag(2,:),'Color',C)
   plot(Vert(1,:),Vert(2,:),'Color',C)
   plot(Hor(1,:),Hor(2,:),'Color',C)
   
   Tx_R = R*[100-i+1;i+1] ;
   Tx_L = R*[-2;100-i+1] ;
   Tx_D = R*[i+1;-3] ;
   
   text(Tx_R(1),Tx_R(2),sprintf('%.0f',i))
   text(Tx_L(1),Tx_L(2),sprintf('%.0f',i),'HorizontalAlignment','right')
   text(Tx_D(1),Tx_D(2),sprintf('%.0f',i))
end

nom = R*[-15 -5 105; 105 -7 -7] ;

text(nom(1,1),nom(2,1),'Metilciclohexano')
text(nom(1,2),nom(2,2),'Heptano')
text(nom(1,3),nom(2,3),'Anilina')
Equi_HC = R*[HC_Anilina';HC_Metilciclohexano'] ;
Equi_S = R*[S_Anilina';S_Metilciclohexano'] ;

for j = 1:13
    plot([Equi_HC(1,j) Equi_S(1,j)],[Equi_HC(2,j) Equi_S(2,j)],'Color',[0.3 0.3 0.3])
end

plot(Equi_HC(1,:),Equi_HC(2,:),'r',Equi_S(1,:),Equi_S(2,:),'b')
axis('off')

%%
figure('Color','White')
tiledlayout(1,2)

nexttile
R = [1 cos(pi/3) ; 0 sin(pi/3)] ;

hold on
for i = 0:10:100
   Diag = [0 i;i 0] ;
   Vert = [i i;0 100-i] ;
   Hor = [0 100-i;i i];
   
   Diag = R*Diag ;
   Vert = R*Vert;
   Hor =  R*Hor;
   
   if i==0 || i==100
   C = [0 0 0] ;
   else
   C = [0.50 0.50 0.50] ;
   end
   
   plot(Diag(1,:),Diag(2,:),'Color',C)
   plot(Vert(1,:),Vert(2,:),'Color',C)
   plot(Hor(1,:),Hor(2,:),'Color',C)
   
   Tx_R = R*[100-i+1;i+1] ;
   Tx_L = R*[-2;100-i+1] ;
   Tx_D = R*[i+1;-3] ;
   
   text(Tx_R(1),Tx_R(2),sprintf('%.0f',i))
   text(Tx_L(1),Tx_L(2),sprintf('%.0f',i),'HorizontalAlignment','right')
   text(Tx_D(1),Tx_D(2),sprintf('%.0f',i))
end

nom = R*[-13 -10 110; 105 -7 -7] ;

text(nom(1,1),nom(2,1),'Metilciclohexano')
text(nom(1,2),nom(2,2),'Heptano')
text(nom(1,3),nom(2,3),'Anilina')
Equi_HC = R*[HC_Anilina';HC_Metilciclohexano'] ;
Equi_S = R*[S_Anilina';S_Metilciclohexano'] ;

for j = 1:13
    plot([Equi_HC(1,j) Equi_S(1,j)],[Equi_HC(2,j) Equi_S(2,j)],'Color',[0.3 0.3 0.3])
end

plot(Equi_HC(1,:),Equi_HC(2,:),'r',Equi_S(1,:),Equi_S(2,:),'b')
axis('off')

nexttile

plot([0 100],[0 100],'k',HC_Metilciclohexano,S_Metilciclohexano,'g')
grid minor
