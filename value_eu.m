
function[matr, vec] = value_eu(N,T,sigma,S0,r,K,option)

N = N;
T = T;
sigma = sigma;
S0 = S0;
r = r;
K=K;

% Buscamos parámetros necesarios
Dt = T/N;
%r = ((1+r)^(1/365))-1; % Convirtiendo de anual a diario
%sigma = sigma/(365^0.5); % Convirtiendo de anual a diario
q = 1/2 + (sqrt(Dt)/(2*sigma))*(r - (sigma^2)/2);
u = exp(sigma * sqrt(Dt));

% Generamos matriz de ceros
precios = zeros(N+1);
nodos = zeros(N+1);

% Creamos matriz de precios (árbol)
for i=0:N
   precios(1,i+1) = u^i * S0;
   for j=1:i
      precios(j+1,i+1) = u^(i-2*j) * S0;
   end
end


% Loop para definir valor de put en T y vector de ejercicio: 1 si se
% ejerce, 0 si no. Depende del tipo de opción (1=put, 2=call).
opti = precios(:,size(precios,2));

if option == 1
    for i=1:length(opti)
        if K > opti(i)
            vopti(i) = K - opti(i);
            voptiTrue(i) = 1;
        else
            vopti(i) = 0;
            voptiTrue(i) = 0;
        end
    end
    
elseif option == 2
    for i=1:length(opti)
    if K < opti(i)
        vopti(i) = opti(i) - K;
        voptiTrue(i) = 1;
    else
        vopti(i) = 0;
        voptiTrue(i) = 0;
    end
    end
end

% Matriz de valor recursiva
delta = exp(-r * Dt);

% Generamos de matriz de valor y sustituimos última fila
valor = zeros(N+1);
valor(:,length(valor))=vopti.';

% Loop para determinar matriz de valor
for i=N:-1:1
for j=N:-1:1
    if i>j
        valor(i,j)=0;
    else
   valor(i,j) =delta*((valor(i,j+1)*q)+(valor(i+1,j+1)*(1-q)));
    end
end
end

% Outputs
matr = valor;
vec = voptiTrue.';