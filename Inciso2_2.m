% INCISO 2.2
% ERNESTO R. PAREDES PÉREZ
close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Este ejercicio consiste en la simulación de proceso AR(1) para ciertos
% parámetros \phi y con un horizonte definido de realizaciones.
% Adicionalmente, se construye el estadístico Dickey-Fuller, para
% comprobación de raíces unitarias.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inciso a) Generando simulación de series

% Plantando semilla
rng(0830);

% Declarando parámetros
T = 500;
K = 1;
Phi = [0.99,1];

% Generando matriz de errores 
e = ones(T,1);
e = normrnd(0,1,size(e));
e2 = normrnd(0,1,size(e));

% Generando series ----------
y1 = zeros(T,1);
    for i=2:T
        y1(i) = Phi(1)*y1(i-1) + e(i);
    end
    
 y2 = zeros(T,1);
    for i=2:T
        y2(i) = Phi(2)*y2(i-1) + e2(i);
    end
 
% Generando estadísticos D-F ----------

% Generando series de lag
for i=2:T
    lay1(i) = y1(i-1);
    lay2(i) = y2(i-1);
end

% Regresiones
lay1 = lay1.';
lay2 = lay2.';

B1 = (lay1'*lay1)\(lay1'*y1);
B2 = (lay2'*lay2)\(lay2'*y2);

% Generando errores de las series
er1 = y1 - y1*B1;
er2 = y2 - y2*B2;

% Encontrando desv estándar de errores
sigma1 = (er1.'*er1)/(T-K);
sigma2 = (er2.'*er2)/(T-K);

% Encontrando varianza de estimadores
var1=sigma1*inv(lay1'*lay1);
var2=sigma2*inv(lay2'*lay2);

tau1 = (B1-1)/sqrt(var1);
tau2 = (B2-1)/sqrt(var2);

 % Gráficos --------------
 plot(y1)
 hold on
 plot(y2)
xlabel('T')
title('Series generadas')
legend('0.99', '1')
legend('Location','northeastoutside')
hold off
saveas(gcf,'II2Series.png')
close

%% Inciso b) Creando función
% Eliminando semilla inicial
rng('default');

% Utilizando función para encontrar parámetros
[phi, tau] = DF(500,1);

%% Inciso c) Realizando subplot

% Declarando parámetros iniciales para loop
m=10^5;
Tm = [25, 50, 100, 250, 500, 1000]

% Generando matrices de almacenamiento
betas=zeros(m, length(Tm));
taus=zeros(m, length(Tm));


% Definimos loops ------

time1= datetime('now');
for i=1:m
    for j=1:length(Tm)
        [phi, tau] = DF(Tm(j), 1);
        betas(i,j)=phi;
        taus(i,j)=tau;
        fprintf('repetición = %f, muestra = %f\n',i,Tm(j))  % Permite que vaya reportando iteración
    end
end
time2=datetime('now');
Duracion=time2 - time1 % Permite saber duración exacta de iteración

% Subplots ------

lfg={'T=25', 'T=50', 'T=100', 'T=250', 'T=500', 'T=1000'};

% Betas --
for i=1:length(Tm)
   subplot(2,3,i) 
    histogram(betas(:,i),'normalization','probability')
    title(lfg(i), 'interpreter', 'latex')
    xlim([-0.45, 1.2]);
    grid on;
end
sgtitle('Histogramas de \beta') 
saveas(gcf,'II2Betas.png')

% Taus --
for i=1:length(Tm)
   subplot(2,3,i) 
    histogram(taus(:,i), 'normalization','probability')
    title(lfg(i), 'interpreter', 'latex')
    grid on;
end
sgtitle('Histogramas de \tau')
saveas(gcf,'II2Taus.png')
close


%% Inciso d) Valores Críticos
percen = [1, 2.5, 5, 10];
vcriticos = zeros(length(percen),size(taus,2));

for i=1:length(percen)
    vcriticos(i,:) = prctile(taus,percen(i));
end

% Exportamos a Latex
rowLabels = {'1\%', '2.5\%', '5\%', '10\%'};
matrix2latex(vcriticos, 'II2.tex', 'rowLabels', rowLabels, 'columnLabels', lfg, 'format', '%0.2f');

%% Inciso e) Evaluación de poder estadístico

% Definimos parámetros
Phis = [0.9, 0.95, 0.98, 0.99]; 
m2=10^4

% Generamos Loop -----

% Por cada uno de los Phis genero una matriz que
for i=1:length(Phis)
  % Por cada repetición de muestra
  for k=1:m2
    % Por cada tamaño de muestra (fila)
   for j=1:length(Tm) 
       % Encuentro el Tau correspondiente dado tamaño (columna)
       [phi, tau] = DF(Tm(j), Phis(i));
       % Si tau obtenido es menor que el valor crítico estimado antes
       if tau < vcriticos(3,j);
           % Se rechaza, pone en matriz
           tausp(k,j) = 1;
       else
           % No se rechaza, pone en matriz
           tausp(k,j) = 0;
       end
       % Comando para dar seguimiento
       fprintf('phi = %f, repetición = %f, muestra = %f\n',Phis(i),k,Tm(j))
   end
  end
  
  % Dado Phi, por cada tamaño de muestra
  for l=1:length(Tm)
      % Pongo en power, porcentaje de rechazos a total de repeticiones
      power(i, l) = (sum(tausp(:,l))/m2)*100;
  end
end



% Exportamos power
rowLabels = {'0.9', '0.95', '0.98', '0.99'};
matrix2latex(power, 'II2Power.tex', 'rowLabels', rowLabels, 'columnLabels', lfg, 'format', '%0.2f');

plot(Tm, power)
legend(rowLabels)
xticks([Tm])
saveas(gcf,'II2Criticos.png')
close

