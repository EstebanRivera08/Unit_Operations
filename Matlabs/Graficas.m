clc
clear

HC_M = [0 9.2 18.6 22 33.8 40.9 46 59 66.8 71.4 ...
    73.3 83.3 88.1]'/100 ;
HC_H = [92.9 83.1 73.4 69.8 57.6 50.4 45 30.7 22.8 ...
    17.8 15.7 5.4 0]'/100 ;
HC_A = 1-HC_M-HC_H ;

S_M = [0 0.8 2.7 3 4.6 6 7.4 9.2 11.3 12.7 13.1 ...
    15.6 16.9]'/100;
S_H = [6.2 6 5.3 5.1 4.5 4 3.6 2.8 2.1 1.6 1.4 ...
    0.6 0]'/100;
S_A = 1-S_M-S_H ;


% CARTESIANAS & EQUILATERO
figure('Color','white')
tiledlayout(1,2)

nexttile
R = [1 0 ; 0 1] ;

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

for j = 1:13
    plot([HC_A(j) S_A(j)],[HC_M(j) S_M(j)],'Color',[0.3 0.3 0.3])
end

nom = [-0.10 -0.10 1; 1.05 -0.07 -0.07] ;

text(nom(1,1),nom(2,1),'Metilciclohexano')
text(nom(1,2),nom(2,2),'Heptano')
text(nom(1,3),nom(2,3),'Anilina')
plot(HC_A,HC_M,'r',S_A,S_M,'b')
axis('off')

nexttile
R = [1 cos(pi/3) ; 0 sin(pi/3)] ;

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

% CARTESIANAS & X,Y
figure('Color','White')
tiledlayout(1,2)

nexttile
R = [1 0 ; 0 1] ;

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

for j = 1:13
    plot([HC_A(j) S_A(j)],[HC_M(j) S_M(j)],'Color',[0.3 0.3 0.3])
end

nom = [-0.10 -0.10 1; 1.05 -0.07 -0.07] ;

text(nom(1,1),nom(2,1),'Metilciclohexano')
text(nom(1,2),nom(2,2),'Heptano')
text(nom(1,3),nom(2,3),'Anilina')
plot(HC_A,HC_M,'r',S_A,S_M,'b')
axis('off')

nexttile 
plot([0 1],[0 1],'k',HC_M,S_M,'g')
xlabel('%wt. Metilciclohexano en Anilina')
ylabel('%wt., Metilciclohexano en Heptano')
grid minor

% Base libre

N = HC_A./(1-HC_A) ;
M = S_A./(1-S_A) ;

X = HC_M./(1-HC_A) ;
Y = S_M./(1-S_A) ;

figure('Color','White')
tiledlayout(2,1)

nexttile

plot(X,N,'r',Y,M,'b')
hold on
for j = 1:13
    plot([X(j) Y(j)],[N(j) M(j)],'Color',[0.3 0.3 0.3])
end
grid minor
xlabel('X,Y')
ylabel('N,M')
xlim([0 1])

nexttile
plot([0 1],[0 1],'k',X,Y,'g')
xlabel('X')
ylabel('Y')
grid minor


