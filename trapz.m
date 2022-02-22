
 function error = trapz(min, max, int, val, nodos)
 
 for n=1:1:nodos               % Número de pasos
    ax=min;                   % Límite inferior de integral
    bx=max;                    % Límite superior de integral
    Valorx = val;      % Valor de ingregral
    h = (bx-ax)/(n-1);          % Tamaño del paso
    eval = ax:h:bx;             % Vector de nodos equidistantes
    integral = int;            % Definimos integral (no pude hacerlo loop `x' como en STATA
    f = integral(eval);      % Evaluamos integral con nodos
 
   
    if n==1
        I_trap1 = h/2*(f(n));
        err(n) = abs(Valorx-I_trap1)/Valorx;
    else
        I_trap1 = h/2*(f(1)+2*sum(f(2:n-1))+f(n));
        err(n) = abs(Valorx-I_trap1)/Valorx;
    end
    
    error=err;
    
    

 end
 end