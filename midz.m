
 function error = midz(min, max, int, val, nodos)
 
 for n=1:1:nodos               % Número de pasos
    ax= min;                   % Límite inferior de integral
    bx=max;                    % Límite superior de integral
    Valorx = val;      % Valor de ingregral
    h = (bx-ax)/(n-1);          % Tamaño del paso
    eval = ax:h:bx;             % Vector de nodos equidistantes
    integral = int;            % Definimos integral (no pude hacerlo loop `x' como en STATA)
    f = integral(eval);      % Evaluamos integral con nodos
   
   % Forma Normal
   
   I_Norm2 = h*sum(f(1:n));
   err(n) = abs(Valorx-I_Norm2)/Valorx;
   
   error=err;
 end
 
 end