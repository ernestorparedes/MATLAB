% INCISO 1
% ERNESTO R. PAREDES PÉREZ
close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% En este inciso procuramos aproximar integrales a partir de tres reglas.
% Así, en la primer sección se programan las integrales. Se intentó
% infructuosamente usar un loop, pero nunca logré encontrar la forma de que
% me guardara las variables de forma análoga al "comando" `x' en STATA para
% los loops. Por eso, se definen vectores que incluyen los valores de ambas
% integrales.
% INCISOS INCLUIDOS (REALIZADOS):
% A) Sección 1
% B) Sección 2 y 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definimos integrales 1 y 2
int1 = @(z) z.^3 - z.^2 - z + 1 ;
int2 = @(z) sin(abs(pi*z));

% Definimos valor real de integrales 1 y 2
TrueVal = [32/3, 4/pi];

% Definimos límites de integrales 1 y 2
a = [-1, -1];
b = [3, 1];

% Definimos rango de pasos
inicio = 1;
fin = 11;

% Matriz para valores almacenados de aproximación

 %% Inciso a) Aproximación de integrales

% Integral 1 -----------------------------------------------------------
v=1;
    
    % Loop para realización de operación
    for n=inicio:1:fin               % Número de pasos
    ax= a(v);                   % Límite inferior de integral
    bx=b(v);                    % Límite superior de integral
    Valorx = TrueVal(v);      % Valor de ingregral
    h = (bx-ax)/(n-1);          % Tamaño del paso
    eval = ax:h:bx;             % Vector de nodos equidistantes
    integral = int1;            % Definimos integral (no pude hacerlo loop `x' como en STATA
    f = integral(eval);      % Evaluamos integral con nodos
   % Iinterval = zeros(n,1);
    
   % Forma Normal -------------
   I_Norm1 = h*sum(f(1:n));
   I1Q1(n) = abs(Valorx-I_Norm1)/Valorx;
   
   if n==fin
   Valor(1,1)=I_Norm1;
   end
   
   % Trapezoidal -------------
    if n==1
        I_trap1 = h/2*(f(n));
        I1Q2(n) = abs(Valorx-I_trap1)/Valorx;
    else
        I_trap1 = h/2*(f(1)+2*sum(f(2:n-1))+f(n));
        I1Q2(n) = abs(Valorx-I_trap1)/Valorx;
        
    if n==fin
    Valor(2,1)=I_trap1;
    end
    end
    

    
   % Simpson -------------
   % Simpson se plantea de esta forma porque me estuvo dando problemas al
   % aproximar. Fui definiendo los valores iniciales para asegurarme que
   % funcionara de una mejor forma, aunque aún mantiene cierta oscilación
   % en los resultados //
    if n==1 % Condición inicial
       I_simp1 = h/3*(f(n));
       I1Q3(n) = abs(Valorx-I_simp1)/Valorx;
    elseif n==2 % Condición pares
       I_simp1 = h/3*(f(1)+f(n));
       I1Q3(n) = abs(Valorx-I_simp1)/Valorx;
    elseif n==3 % Condición impares
        I_simp1 = h/3*(f(1)+4*sum(f(2:2:n-1))+f(n));
       I1Q3(n) = abs(Valorx-I_simp1)/Valorx;
    else
    I_simp1 = h/3*(f(1)+4*sum(f(2:2:n-1))+2*sum(f(3:2:n-1))+f(n));
       I1Q3(n) = abs(Valorx-I_simp1)/Valorx;
    end
    
   if n==fin
   Valor(3,1)=I_simp1;
   end
    end

% Integral 2 -----------------------------------------------------------
v=2; 

    % Loop para realización de operación
    for n=inicio:1:fin               % Número de pasos
    ax= a(v);                   % Límite inferior de integral
    bx=b(v);                    % Límite superior de integral
    Valorx = TrueVal(v);      % Valor de ingregral
    h = (bx-ax)/(n-1);          % Tamaño del paso
    eval = ax:h:bx;             % Vector de nodos equidistantes
    integral = int2;            % Definimos integral (no pude hacerlo loop `x' como en STATA
    f = integral(eval);      % Evaluamos integral con nodos
   
   % Forma Normal -------------
   I_Norm2 = h*sum(f(1:n));
   I2Q1(n) = abs(Valorx-I_Norm2)/Valorx;
       if n==fin
       Valor(1,2)=I_Norm2;
       end
   % Trapezoidal -------------
    if n==1
        I_trap2 = h/2*(f(n));
        I2Q2(n) = abs(Valorx-I_trap2)/Valorx;
    else
        I_trap2 = h/2*(f(1)+2*sum(f(2:n-1))+f(n));
        I2Q2(n) = abs(Valorx-I_trap2)/Valorx;
       if n==fin
       Valor(2,2)=I_trap2;
       end
    end
    
   % Simpson -------------
    if n==1
       I_simp2 = h/3*(f(n));
       I2Q3(n) = abs(Valorx-I_simp2)/Valorx;
    elseif n==2
       I_simp2 = h/3*(f(1)+f(n));
       I2Q3(n) = abs(Valorx-I_simp2)/Valorx;
    elseif n==3
        I_simp2 = h/3*(f(1)+4*sum(f(2:2:n-1))+f(n));
       I2Q3(n) = abs(Valorx-I_simp2)/Valorx;
    else
    I_simp2 = h/3*(f(1)+4*sum(f(2:2:n-1))+2*sum(f(3:2:n-1))+f(n));
       I2Q3(n) = abs(Valorx-I_simp2)/Valorx;
    end
       if n==fin
       Valor(3,2)=I_simp2;
       end
    end 
    
 % Exportando Tablas
 rowLabels = {'Midz', 'Trapz', 'Simps'};
 columnLabels = {'1', '2'};
 I1 = [I1Q1; I1Q2; I1Q3];
 I2 = [I2Q1; I2Q2; I2Q3];
 matrix2latex(I1, 'I1.tex', 'rowLabels', rowLabels, 'format', '%0.2f');
 matrix2latex(I2, 'I2.tex', 'rowLabels', rowLabels, 'format', '%0.2f'); 
 matrix2latex(Valor, 'IV.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'format', '%0.2f'); 

 
 %% Inciso b) Creación de funciones
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % En este inciso convertimos las reglas utilizadas antes en funciones
 % individuales midz, trapz y simps, para cada regla correspondiente. Se
 % toman en cuenta 11 nodos.
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 nds=11;
 
 I1_midz = midz(-1,3,int1,32/3,nds);
 I1_trapz = trapz(-1,3,int1,32/3,nds);
 I1_simps = simps(-1,3,int1,32/3,nds);
 
 I2_midz = midz(-1,3,int2,4/pi,nds);
 I2_trapz = trapz(-1,3,int2,4/pi,nds);
 I2_simps = simps(-1,3,int2,4/pi,nds);

 
%% Gráficos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% En esta sección procedemos a realizar los gráficos.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Integral 1
nexttile
    plot(I1_midz)
    xlabel('Nodos')
    ylabel('Error')
    title('Integral 1')
hold on
    plot(I1_trapz)
    plot(I1_simps)
hold off
    legend('Midz', 'Trapz', 'Simps')

% Integral 2
nexttile
    plot(I2_midz)
    xlabel('Nodos')
    ylabel('Error')
    title('Integral 2')
hold on
    plot(I2_trapz)
    plot(I2_simps)
hold off
    legend('Midz', 'Trapz', 'Simps')
saveas(gcf,'IA.png')
close
