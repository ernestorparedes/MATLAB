% INCISO 2.1
% ERNESTO R. PAREDES PÉREZ
close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Este inciso consiste en la extracción de ciclo y tendencias para series,
% con tres métodos. A saber: log-lineal, log-cuadrático y Hodrick-Prescott.
% Particularmente para *Hodrick-Prescott es necesario hacer una aclaración:
% Al buscar la definición de matriz tridiagonal tomé la idea que tenía que
% ser cuadrada y con las diagonalas ordenas como 1, -2, 1. Sin embargo,
% para que funcionara correctamente, fue necesario realizarle un par de
% cambios a dicha matriz en los 4 valores iniciales de la diagonal y los
% 4 últimos. Los resultados son similares a utilizar la especificación de
% la tarea (matriz n x (n-2)).
% INCISOS REALIZADOS: TODOS.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Al correrlo da un error, pero no tiene efectos prácticos. El problema
% está con la regresión, sin embargo usé el comando regress y da resultados
% iguales.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Inciso a) Importando base de datos

filename = 'BaseDatos.csv';
delimiterIn = ',';
headercolIn = 1;
y = importdata(filename, delimiterIn, headercolIn );
y=y.data;
years = y(:,1);
y = log(y);
y(:,[1]) = [];

Chile = y(:,1);
Chile = reshape(Chile,[],1);
Arg = y(:,2);

%% Inciso b) Extracción de ciclo y tendencia

% Hodrick-Prescott ---------------------------------

% Definimos tamaño de serie
n = length(Chile);

% Generando matriz Identidad
I = eye(n);

% Generando matriz K
a = -2; % Columna de en medio
b = 1; % Columna de abajo
c = 1; % Columna de arriba
K = diag(a*ones(1,n)) + diag(b*ones(1,n-1),1) + diag(c*ones(1,n-1),-1);

% Generando matriz multiplicada
KK = K.'*K;

% Aquí se programa la matriz descrita en la tarea
K2 = diag(b*ones(1,n)) + diag(a*ones(1,n-1),1) + diag(c*ones(1,n-2),2);
K2 = K2(1:length(K)-2,:);
KK2 = K2.'*K2;

% Cambiando valores para mantener homogeneidad
KK(1,1)=6;
KK(n,n)=6;

% Definimos vector Lambda (contendrá todos los valores)
Lb = [0, 6.25, 10, 20, 100, 225, 500, 750, 900];

% Chile -----------
 
% Loop para generación de series con cada valor de Lb
for i = 1:length(Lb)

% Generamos Matriz A
A = I + Lb(i)*KK ;
% Realizamos cambios en ciertos valores
A(1,1)=1+Lb(i);       A(1,2)=-2*Lb(i);
A(2,1)=-2*Lb(i);      A(2,2)=5*Lb(i)+1;
A(n-1,n-1)=5*Lb(i)+1; A(n-1,n)=-2*Lb(i);
A(n,n-1)=-2*Lb(i);    A(n,n)=1+Lb(i);

% En este punto KK y KK2 son iguales

% Matriz que contiene tendencia
Chs = inv(A)*Chile;
ChileS(i,:) = Chs';

% Matriz que contiene ciclo
Chc = (I-inv(A))*Chile;
ChileC(i,:) = Chc';
end

% Tratamos matrices y graficamos
ChileC = ChileC.';
ChileS = ChileS.';
years=years.';


% Argentina -----------

% Loop para generación de series con cada valor de Lb
for i = 1:length(Lb)

% Generamos Matriz A
A = I + Lb(i)*KK ;

% Realizamos cambios en ciertos valores
A(1,1)=1+Lb(i);       A(1,2)=-2*Lb(i);
A(2,1)=-2*Lb(i);      A(2,2)=5*Lb(i)+1;
A(n-1,n-1)=5*Lb(i)+1; A(n-1,n)=-2*Lb(i);
A(n,n-1)=-2*Lb(i);    A(n,n)=1+Lb(i);

% Matriz que contiene tendencia
Args = inv(A)*Arg;
ArgS(i,:) = Args';

% Matriz que contiene ciclo
Argc = (I-inv(A))*Arg;
ArgC(i,:) = Argc';
end

% Tratamos matrices y graficamos
ArgC = ArgC.';
ArgS = ArgS.';
years=years.';

% Log-Linear Detrending ---------------------------------

C = ones(n,1);
Z = [C, years];

%%%% Chile

B = (Z.'*Z)\(Z.'*Chile);
CLLc = Chile - Z*B; % Ciclo
CLLt = Z*B;         % Tendencia

%%%% Arg
B = (Z.'*Z)\(Z.'*Arg);
ALLc = Arg - Z*B;   % Ciclo
ALLt = Z*B;         % Tendencia

% Log-Quadratic Detrending ---------------------------------

C = ones(n,1);
Z = [C, years, years.^2];
%%%% Chile
B = (Z.'*Z)\(Z.'*Chile);
CLQc = Chile - Z*B; % Ciclo
CLQt = Z*B;         % Tendencia

%%%% Arg
B = (Z.'*Z)\(Z.'*Arg);
ALQc = Arg - Z*B;   % Ciclo
ALQt = Z*B;         % Tendencia

%% Inciso c) Gráficos 

CArg = [ALLc, ALQc];
TArg = [Arg, ALLt, ALQt];

CChl = [CLLc, CLQc];
TChl = [Chile, CLLt, CLQt];

% Argentina ----------------

% Ciclo --
nexttile
plot(years, CArg)
xlabel('Años')
ylabel('% desviación de la tendencia')
title('Argentina: Ciclo')
legend('Linear', 'Quadratic')
legend('Location','northeastoutside')

% Tendencia --
nexttile
plot(years, TArg)
xlabel('Años')
ylabel('Log PIB per cápita')
title('Argentina: Tendencia')
legend('Original', 'Linear', 'Quadratic')
legend('Location','northeastoutside')
saveas(gcf,'IICACT.png')
close

% Hodrick-Prescott --
nexttile
plot(years, ArgC')
xlabel('Años')
ylabel('% desviación de la tendencia')
title('Argentina: Ciclo')
legend('Original (=0)', '6.25', '10', '20', '100', '225', '500', '750', '900');
legend('Location','northeastoutside')

nexttile
plot(years, ArgS')
xlabel('Años')
ylabel('Log PIB per cápita')
title('Argentina: Tendencia')
legend('Original (=0)', '6.25', '10', '20', '100', '225', '500', '750', '900');
legend('Location','northeastoutside')
saveas(gcf,'IICAHP.png')
close 
% Chile ----------------

% Ciclo --
nexttile
plot(years, CChl)
xlabel('Años')
ylabel('% desviación de la tendencia')
title('Chile: Ciclo')
legend('Linear', 'Quadratic')
legend('Location','northeastoutside')

% Tendencia -- 
nexttile
plot(years, TChl)
xlabel('Años')
ylabel('% Tendencia')
title('Chile: Tendencia')
legend('Original', 'Linear', 'Quadratic')
legend('Location','northeastoutside')
saveas(gcf,'IICCCT.png')
close 

% Hodrick-Prescott --
nexttile
plot(years, ChileC')
xlabel('Años')
ylabel('% desviación de la tendencia')
title('Chile: Ciclo')
legend('Original (=0)', '6.25', '10', '20', '100', '225', '500', '750', '900');
legend('Location','northeastoutside')

nexttile
plot(years, ChileS')
xlabel('Años')
ylabel('Log PIB PC')
title('Chile: Tendencia')
legend('Original (=0)', '6.25', '10', '20', '100', '225', '500', '750', '900');
legend('Location','northeastoutside')
saveas(gcf,'IICCHP.png')
close