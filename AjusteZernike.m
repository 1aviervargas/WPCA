function [S, Zzer, Zpol] = AjusteZernike(X, Y, Z, Mask, n)


% S = AjusteZernike(R, Theta, Z, n)
% [S, Zzer] = AjusteZernike(R, Theta, Z, n)
%
% Esta funci�n ajusta unos datos de sagita (R, Theta, Z) a n polinomios de
% Zernike, utilizando la funci�n generadora Zernike. Los polinomios son no
% normalizados, y se utiliza la numeraci�n de Wyant.
%
% Los argumentos de salida son S (vector de coeficientes de longitud n+1) y
% opcionalmente Zzer, que son los valores de sagita estimados mediante la 
% expansi�n polinomial en los puntos de coordenadas (R, Theta)
%
% Al utilizarse un muestreo de los polinomios y funciones que se expanden,
% los coeficientes no pueden obtenerse directamente por proyecci�n. En su
% lugar se resuelve el problema en t�rminos de m�nimos cudradados,
% utilizando la inversa generalizada de Penrose:
% inv(A'*A)*A'

% Autor: Jos� Alonso
% Versiones:    1.0 Argumentos de salida: coeficientes y polinomios
%               2.0 Argumentos de salida: coeficientes y, opcionalmente, la
%               estimaci�n de la superficie con la expansi�n.
% Fecha: Julio de 2005
% _________________________________________________________________________


    maxX = max(X(:));
    maxY = max(Y(:));
    
    mX = mean2(X);
    mY = mean2(Y);
    
    Xn = (X - mX)./maxX;
    Yn = (Y - mY)./maxY;
    
    [Theta,R]=cart2pol(Xn,Yn);
    R = R(:);
    Theta = Theta(:);
    
    % Generaci�n de polinomios
    R=R./(max(R)+eps);
    Zpol = Zernike(n, R, Theta);


% C�lculo de coeficientes por m�nimos cuadrados (inversa generalizada de
% Penrose)
S = inv(Zpol(Mask(:),:)'*Zpol(Mask(:),:))*Zpol(Mask(:),:)'*Z(Mask(:));

% C�lculo de la superficie estimada
if nargout == 3
    Zzer = Zpol*S;
end

Zzer = reshape(Zzer,size(Z,1),size(Z,2));