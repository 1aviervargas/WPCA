function [varargout] = pcaDemod(varargin)
%PCADEMOD This function obtains the wrapped phase from a sequence of
%phase-shifted interferograms using the Principal Component Analysis
%algorithm (PCA). This function is inspired in the work of A. Leonardis,
%D. Skocaj that is available from CMP Vision Algorithms
%http://vicos.fri.uni-lj.si/danijels/downloads
%http://visionbook.felk.cvut.cz
%
% PCA is a linear integral transformation that simplifies
% a multidimensional dataset to a lower dimension.
% The implementation of the function pca
% uses the efficient implementation of singular
% value decomposition (svd).
%
% Usage: [pw,Mod,U1,U2,V] = pca(I,K,Mask)
% Inputs:
%   I  [NRows x NCols x num] Image 3D matrix where NRows and NCols are the 
%   number of rows and columns of the fringe patterns and num is the number 
%   of interferograms, I(:,:,1) is the first interfergramA.   
%   K  1x1  Number of eigenvalues and eigenvectors to be returned.
%   Mask [NRows x NCols] is the processing mask
% Outputs:
%   pw  [NRows x NCols]  Wrapped modulating phase.
%   U1  [NRows x NCols]  First principal component eigenvector.
%   U2  [NRows x NCols]  Second principal component eigenvector.
%   V  [NRows x NCols]  Vector with the differents eigenvalues.
%   Javier Vargas, 
%   25/11/10 
%   Copyright 2010 
%   LINES-INTA 
%   $ Revision: 1.0.0.0 $
%   $ Date: 25/10/10 $

try
    funcName='pcaDemod';
    numvarargs = length(varargin);
    if numvarargs > 3
        error(['MyToolbox:' funcName ':TooManyInputs'], ...
            'requires at most 2 optional inputs: compound interferogram data and number of eigenvalues');
    end
    
    dim = size(varargin{1});
    %Compound image formed by the different interferograms columnwise
    %stacked
    X=reshape(varargin{1},dim(1)*dim(2),dim(3));
    %Number of eigenvalues and eigenvectors to be returned.
    K = varargin{2};
    Mask = varargin{3};
    
    for i=1:dim(3)
        temp = X(:,i);
        %Compound image after filtering the mask
        Xf(:,i)= temp(Mask(:));
    end
    
    
    [M N]=size(Xf);
    Xm=mean(Xf,2);
    Xd=Xf-repmat(Xm,1,N);
    if (N < M) %less images than image length
        C=Xd'*Xd;
        [V D Vt]=svd(C);
        U=Xd*V;
        U=U./repmat(sqrt(diag(D)'),M,1);
    else %more images than image length
        C=Xd*Xd';
        [U D Ut]=svd(C);
    end;
    
    L=diag(D)'/N;    
    if nargin>1
        U=U(:,1:K);
%        L=L(1:K);
    end;    

    for i=2:K            
        U(:,i)=max(U(:,1)).*(U(:,i)./max(U(:,i)));
    end
    
    %U(:,1)=max(U(:,2)).*(U(:,1)./max(U(:,1)));
    %U1 = X(:,1).*0;
    %U2 = U1;
    T = double(Mask);
    for i=1:K
        T(Mask) = U(:,i);
        PC(:,:,i) = T;
    end
    
    pw = reshape(atan2(PC(:,:,1),PC(:,:,2)),dim(1),dim(2));
    Mod =  reshape(sqrt(PC(:,:,1).^2+PC(:,:,2).^2),dim(1),dim(2));
    
    T(Mask) = Xm;
    
    varargout{1} = pw;
    varargout{2} = Mod;
    varargout{3} = PC;
    varargout{4} = L;
    varargout{5} = T;
    
catch ME
    throw(ME)
end
