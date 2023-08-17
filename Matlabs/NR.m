function f=NR(f,df,p0)

    it = 0;                     % Creamos una variable para controlar un máximo de iteraciones
 
    
    while abs(f(p0))>1e-10 && it<10000      % Se hace mientras el valor absoluto de f(p0)(posible raiz) sea mayor al error
        p = p0-(f(p0)/df(p0));      % Calculamos el siguiente punto, el cual converge más cerca a la raiz
   
        p0 = p;                 % Nuestro siguiente punto de referencia es el punto que acabamos de calcular
        it = it + 1;            % Actualizamos el número de iteraciones
    end
    
    f = p0;
    