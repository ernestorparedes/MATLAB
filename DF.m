

function[phi, tau] = DF(T, ph)

%T=100;
e = ones(T,1);
e = normrnd(0,1,size(e));

% Generando serie
y1 = zeros(T,1);
    for i=2:T
        y1(i) = ph*y1(i-1) + e(i);
    end
    
% Generando serie de lag
for i=2:T
    lay1(i) = y1(i-1);
end

% Regresión
lay1 = lay1.';
B1 = (lay1'*lay1)\(lay1'*y1);

% Generando errores y desviación estándar de estos
er1 = y1 - B1*lay1;
sigma1 = (er1.'*er1)/(T-1);

% Encontrando varianza de estimador
var=sigma1*inv(lay1'*lay1);

% Encontrando estadístico DF
tau1 = (B1-1)/sqrt(var);

% Outputs
phi=B1;
tau=tau1;