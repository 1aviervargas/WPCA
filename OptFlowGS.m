function [U,V] = OptFlowGS(A,B,K,epsilon,lambda)

%OPTFLOW obtains the optical flow between images A and B using the method presented in [1]. 
%The method is based on a regularized optical flow algorithm [1]. 
%
% Usage: [U,V] = rofDemod(A,B,K,epsilon,lambda)
%
% Inputs:
%
%   A  [NRows x NCols] Image map with NRows and NCols the number of rows 
%       and columns of the fringe pattern
%   B  [NRows x NCols] Image map with NRows and NCols the number of rows 
%       and columns of the fringe pattern   
%   K   [1x1 double value] Max number of iterations of the regularization
%       process. (Default ~ 5-10)
%   elpsilon [1x1 double value] Stopping criteria of the regularization
%       process. (Default ~ 0.3)
%   lambda [1x1 double value] Regularization parameter (Default ~ 10)
%
% Outputs:

%   U  [NRows x NCols]  Derivation of the phase with respect to the X coordinate.
%   V [NRows x NCols]  Derivation of the phase with respect to the Y coordinate.
%
%   Javier Vargas, 
%   06/06/2011 
%   copyright @2011 
%   Biocomputing Unit, Centro Nacional de Biotecnología (CSIC)
%   http://biocomp.cnb.csic.es/
%   $ Revision: 1.0.0.0 $
%   $ Date: 25/10/10 $
%   $ Revision: 1.0.0.0 $
%   $ Date: 06/06/11 $
%
% THIS CODE IS GIVEN WITHOUT ANY GUARANTY AND WITHOUT ANY SUPPORT
%
% REFERENCES 

try
    % Filtrado pre-processing para reducir ruido
    sx = 5;
    sy = 5;
    % Fijamos el tamaño de bloque a procesar
    %bx=8; by=8;
    bx=15; by=15;
    
    [M,N,d] = size(A);
    
    im1 = double(A);
    im2 = double(B);
    
    h = ones(sx,sy)/(sx*sy);
    
    A0 = conv2(im1(:,:,1),h,'same');
    A1 = conv2(im2(:,:,1),h,'same');
    
    % Tamaño de imagen
    nf=M; nc=N;
    
    
    % inicializamos las matrices del flujo óptico
    U = zeros(M,N);
    V = U;
    
    % Calculamos las derivadas
    Hx=[-1 0 1]; Hy=[-1 0 1]';
    
    fx=imfilter(A0,Hx); fy=imfilter(A0,Hy);
    ft=(A1-A0);
    
    nucleo = ones(bx,by);
    fx2 = fx.^2; fy2 = fy.^2; ft2 = ft.^2;
    D = lambda^2 + fx2 + fy2;
    k = 1;
    while k <= K
        % estimar el flujo óptico medio en un entorno de vecindad
        u_medio = imfilter(U,nucleo,'same');
        v_medio = imfilter(V,nucleo,'same');
        
        P = fx.*u_medio + fy.*u_medio + ft;
        U = u_medio - fx.*P./D;
        V = v_medio - fy.*P./D;
        
        % obtener la energía total
        Ux = imfilter(U,Hx);
        Ux2 = Ux.^2;
        Uy = imfilter(U,Hy);
        Uy2 = Uy.^2;
        Vx = imfilter(V,Hx);
        Vx2 = Vx.^2;
        Vy = imfilter(V,Hy);
        Vy2 = Vy.^2;
        
        E = (sum(sum((fx.*U + fy.*V + ft).^2 + lambda*(Ux2 + Uy2 + Vx2 + Vy2)))).^2;
        
        c1 = 'Iteración '; c2 = num2str(k);
        %disp([c1,c2]);
        
        %disp('Energía: '); disp(E);
        if E < epsilon
            break
        end
        k = k + 1;
    end % bucle que especifica el máximo número de iteraciones
    
catch Exception
    rethrow(Exception);
end
