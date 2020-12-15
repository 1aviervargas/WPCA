function Z = Zernike(p, R, Theta);

% Z = Zernike(A,X,Y)
%
% Esta funci�n calcula lpolinomios de Zernike de grado arbitrario de
% acuerdo con la numeraci�n de Wyant-Creatch ("Basic wavefront aberration
% theory for optical metrology" in Applied Optics and Optical Engineering,
% vol. XI). Los argumentos de entrada son:
%
% p -> vector de n�meros enteros con los �ndices de los polinomios que
% desean calcularse. Normalmente es un vector del tipo 0:Np, pero no es
% necesario, ya que pueden ser polinomios aislados
%
% R, Theta -> Vectores columna con las coordenadas polares de los puntos
% en los que se desea calcular los polinomios. Normalmente, se obtinen de
% la forma [Theta, R] = cart2polar(X,Y), X = X(:); Y = Y(:);

% Los polinomios de Zernike son ortogonales solo en el c�rculo unidad y con
% un producto escalar continuo. Cuando se utiliza un producto escalar
% discreto, puede ocurrir que los polinomios de Zernike dejen de ser
% ortogonales. La funci�n Zernike puede calcular los polinomios sobre
% cualquier conjunto de puntos, pero si la coordenada R excede la unidad se
% genera un mensaje de aviso.

% El argumento de salida es una matriz Z en la que la columna i-�sima
% contine los valores del polinomio i-�simo calculado en los puntos
% (R, Theta).


% Autor: Jos� Alonso
% Versiones:    1.0 Ordenaci�n de polinomios de Born&Wolf
%               2.0 Ordenaci�n de polinomios de Wyant. C�digo m�s compacto
%               y sin chequeos.
% Fecha: Julio de 2005
% _________________________________________________________________________


Nz = length(p); Np = length(R);


%__________________________________________________________________________
% Comprobaci�n de la validez del rango
I = find(R > 1);
if ~isempty(I)
    warning('Coordenadas fuera del c�rculo unidad')
end
    

%__________________________________________________________________________
% C�lculo de los �ndices n,m y paridad de los polinomios de Zernike
n = floor(sqrt(p));
m = (n - 1) - floor((p - n.*n)/2 -1);
l = (n+1).^2 - p - 1;
paridad = sign(l.*((1-2*l)/4 + floor(l/2)));

%__________________________________________________________________________
% C�lculo de las potencias de R
% C�lculo de exponentes
ExponentesR = [];
for i = 1:Nz
    s = 0:(n(i)-m(i));
    ExponentesR = [ExponentesR, 2*(n(i)-s)-m(i)];
end
ExponentesR = sort(unique(ExponentesR));
% Prec�lculo de las potencias de R
PotenciasR = ones(Np, length(ExponentesR));
PotenciasR(:,1) = R.^ExponentesR(1);
for i = 2:length(ExponentesR)
    dExp = ExponentesR(i)-ExponentesR(i-1);
    if dExp == 1
        PotenciasR(:,i) = PotenciasR(:,i-1).*R;
    else
        PotenciasR(:,i) = PotenciasR(:,i-1).*(R.^dExp);
    end
end

%__________________________________________________________________________
% C�lculo de los polinomios de Zernike
Z = zeros(Np, Nz);
for i = 1:Nz
    nm = n(i) - m(i);
    for s = 0:nm
        PotenciaR = 2*(n(i) - s) - m(i);
        iPot = find(ExponentesR == PotenciaR);
        c = fct(n(i) + nm - s)/(fct(s)*fct(n(i)-s)*fct(nm-s));
        signo = 1 - 2*mod(s,2);
        Z(:,i) = Z(:,i) + signo * c * PotenciasR(:,iPot);
    end
    if paridad(i) == 1
        Z(:,i) = Z(:,i).*cos(m(i)*Theta);
    elseif paridad(i) == -1
        Z(:,i) = Z(:,i).*sin(m(i)*Theta);
    end
end


        
%__________________________________________________________________________        
% Sub-Funci�n Factorial (solo para escalares)

function y = fct(x)

if x == 0
    y = 1;
else
    y = prod(1:x);
end
    
