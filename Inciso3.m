% INCISO 3
% ERNESTO R. PAREDES PÉREZ
close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% En este inciso se realizan valuación de opciones put y call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definiendo parámetros del modelo
% Definiendo parámetros
N = 180;
T = 0.5;
sigma = 0.2;
S0 = 2;
r = 0.05;
K=2.1;

% Buscamos parámetros u y q
Dt = T/N;
%r = ((1+r)^(1/365))-1; % Convirtiendo de anual a diario
%sigma = sigma/(365^0.5); % Convirtiendo de anual a diario
q = 1/2 + (sqrt(Dt)/(2*sigma))*(r - (sigma^2)/2);
u = exp(sigma * sqrt(Dt));


%% Inciso a) Creación de árbol de precios

% Generamos matriz de ceros para precios y nodos
precios = zeros(N+1);
probs = zeros(N+1);
nodos = zeros(N+1);

% Loop para matriz de precios (árbol)
for i=0:N
   precios(1,i+1) = u^i * S0;
   for j=1:i
      precios(j+1,i+1) = u^(i-2*j) * S0;
   end
end

% Loop para matriz de probabilidades acumuladas (No se utiliza, pero al
% inicio pensé que sí jeje)
for i=0:N
   probs(1,i+1) = q;
   for j=1:i
      if j==i
      probs(j+1,i+1) = 1-q;
      else
      probs(j+1,i+1) = (q)*(1-q);
      end
   end
end

%% Inciso b) Vector de ejercicio

% Considenrando una put, entonces V_t = {max(K-St), 0}
% Seleccionamos únicamente el último vector

puts = precios(:,size(probs,2));

% Loop para definir valor de put en T y vector de ejercicio: 1 si se
% ejerce, 0 si no
for i=1:length(puts)
    if K > puts(i)
        vputs(i) = K - puts(i);
        vputsTrue(i) = 1;
    else
        vputs(i) = 0;
        vputsTrue(i) = 0;
    end
end

%% Inciso c) Matriz de valor recursiva

% Tenemos que construir la matriz con base en la fórmula de valor
delta = exp(-r * Dt);

% Generamos de matriz de valor y sustituimos última fila
valor = zeros(N+1);
valor(:,length(valor))=vputs.';

% Loop para determinar matriz de valor
for i=N:-1:1
for j=N:-1:1
    if i>j
        valor(i,j)=0;
    else
   valor(i,j) =delta*((valor(i,j+1)*q)+(valor(i+1,j+1)*(1-q)));
   %valor(i,j)
   %=delta*((valor(i,j+1)*(probs(i,j+1)))+(valor(i+1,j+1)*(probs(i+1,j+1))));
   %(Esto no se ocupa)
    end
end
end

%% Inciso d) Función value_ue
% Incluye los valores N,T,sigma,S0,r,K y option, respectivamente. Option es
% 1 si put y 2 si call.

[pmatriz, pvec] = value_eu(180,0.5,0.2,2,0.05,2.1,1);

%% Inciso e) Recorridos

r=5000;

% Sembramos matriz y definimos matriz de recorridos
rng(1979);
% La matriz contiene 1 y 0, correspondiente a subida o bajada
recorridos = randi([0 1],r,N);
rprecios = zeros(r,N+1);
% El precio inicial corresponde S0
rprecios(:,1) = S0;

for i=1:r
    for j=2:size(recorridos,2)
        % Si correspondió 1, entonces se multiplica por u
        if recorridos(i,j) == 1
            rprecios(i,j) = rprecios(i,j-1)*u;
        else
         % Si correspondió 0, entonces se divide por u
            rprecios(i,j) = rprecios(i,j-1)/u;
        end
    end
end

% Tomamos precios en T
vcall = rprecios(:,N);


% Generamos loop en vector precios T para definir valor y ejercicio
for i=1:length(vcall)
    % Definimos regla de ejercicio y almacenamiento
    if vcall(i) > K
        vcallv(i) = vcall(i) - K;
        % Si se ejerce, toma valor 1
        vcalle(i) = 1;
    else
        vcallv(i) = 0;
        % Si no se ejerce, toma valor 0
        vcalle(i) = 0;
    end
end

% Calculamos media y proporción
vcallm = mean(vcallv);
vcalls = sum(vcalle);
vcallp = vcalls/r*100;

% Realizamos histograma
nexttile
 histogram(vcallv)
nexttile
 histogram(vcallv, 'normalization','probability')
 sgtitle('Histogramas ejecución de opción')
saveas(gcf,'IIIHis.png')
close
