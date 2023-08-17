   
function PuntoFalso = PF(f,a,b)

    if f(a)*f(b)>0
        error('En el intervalo ingresado no hay una única raíz, ingrese otro')
    else
            c=b-f(b)*(b-a)/(f(b)-f(a)) ; %punto inicial
    k = 0 ;
        while (abs(f(c))>1e-10 && k<1000)    

               if (f(a)*f(c)>0)        
                 a=c ;
                end 

                if (f(a)*f(c)<0)        
                 b=c ;        
                end

           c=b-f(b)*(b-a)/(f(b)-f(a)) ; % Redefinir C
           k=k+1 ;    %#iteración
        
        end
    end
    f
    c
    f(c)
    PuntoFalso = c ;
    
    end  
