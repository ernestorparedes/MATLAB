

 function error = simps(min,max,int,val,nodos)
 
 for n=1:1:nodos               % Número de pasos
    ax=min;                   % Límite inferior de integral
    bx=max;                    % Límite superior de integral
    Valorx = val;      % Valor de ingregral
    h = (bx-ax)/(n-1);          % Tamaño del paso
    eval = ax:h:bx;             % Vector de nodos equidistantes
    integral = int;            % Definimos integral (no pude hacerlo loop `x' como en STATA
    f = integral(eval);      % Evaluamos integral con nodos
 
   % simp = h/3*(f(1)+ 2*sum(f(3:2:n-2))+4*sum(f(2:2:n-1))+f(n));
   % err(n) = abs(Valorx-simp);
   % error = err;
    
   
    if n==1
       simp = h/3*(f(n));
       err(n) = abs(Valorx-simp)/Valorx;
    elseif n==2
       simp = h/3*(f(1)+f(n));
       err(n) = abs(Valorx-simp)/Valorx;
    elseif n==3
        simp = h/3*(f(1)+4*sum(f(2:2:n-1))+f(n));
       err(n) = abs(Valorx-simp)/Valorx;
    else
    simp = h/3*(f(1)+4*sum(f(2:2:n-1))+2*sum(f(3:2:n-1))+f(n));
       err(n) = abs(Valorx-simp)/Valorx;
    end
   
    error = err;
 end
 end