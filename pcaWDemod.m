function [varargout] = pcaWDemod(varargin)
%PCAWDEMOD This function obtains the wrapped phase from a sequence of
%phase-shifted interferograms using the Weighted Principal Component Analysis
%algorithm (PCA).
%
% Usage: [pw,Mod,U1,U2,V] = pcaWDemod(I,S,Mask)
% Inputs:
%   I  [NRows x NCols x num] Image 3D matrix where NRows and NCols are the 
%   number of rows and columns of the fringe patterns and num is the number 
%   of interferograms, I(:,:,1) is the first interfergramA.   
%   K  1x1  Number of eigenvalues and eigenvectors to be returned.
%   Mask [NRows x NCols] is the processing mask
%   Seed [NRows x NCols x 2] First guess of the 2 Principal Components
%   obtained by a two-step demodulation approach.
%   Alpha  1x1  Significance value. Gives more weight to good interferograms
%   when is close to 0 and less weight when is close to 1. Values between
%   [0.5-1] should work well.
% Outputs:
%   pw  [NRows x NCols]  Wrapped modulating phase.
%   Mod [NRows x NCols]  Modulation term phase.
%   PC  [NRows x NCols x K]  Principal component eigenvectors.
%   L   [1 x K]  Eigenvalues.
%   w   [1 x num] weights used in the wPCA approach
%   Javier Vargas
%   25/06/20 
%   Copyright 2020
%   Universidad Complutense de Madrid 
%   $ Revision: 1.0.0.0 $
%   $ Date: 25/06/20 $

try
    funcName='pcaWDemod';
    numvarargs = length(varargin);
    if numvarargs > 5
        error(['MyToolbox:' funcName ':TooManyInputs'], ...
            'requires at most 5 inputs');
    end        

    %Interferograms
    I = varargin{1};
    %Number of eigenvectors/eigenvalues to compute
    K = varargin{2};
    %Mask
    Mask = varargin{3};
    %Initial seed
    seed = varargin{4};
    %Initial seed
    alpha = varargin{5};
    
    if (ischar(seed))
        if (strcmp(seed,'pca'))
            [pw,Mod,PC]=pcaDemod(I,K,Mask);
        else    
            [pw,Mod,PC] = gsDemod(I(:,:,1),I(:,:,2),Mask);
        end
    else        
        PC = seed;
        pw = atan2(-PC(:,:,1),PC(:,:,2));
    end
           
    Xm = sum(I,3)/size(I,3);
    
    %Compound image formed by the different interferograms columnwise
    %stacked
    dim = size(varargin{1});
    X=reshape(varargin{1},dim(1)*dim(2),dim(3));
    
    for i=1:dim(3)
        temp = X(:,i);
        %Compound image after filtering the mask
        Xf(:,i)= temp(Mask(:));
    end
    
    [M N]=size(Xf);
    
    for i= 1:2
        
        A = is2fs( I, Xm, PC, 2, Mask);        
        w   =  (A(1,:).^2+A(2,:).^2);
        w   =  (w./max(w(:))).^alpha;
        
        [U,A,D,tsquared,explained,mu] = pca(Xf','weight',w);
        U=U./repmat(sqrt((D)'),M,1);
        
        L=diag(D)'/N;
        if nargin>1
            U=U(:,1:K);
        end;
        
        for i=2:K
            U(:,i)=max(U(:,1)).*(U(:,i)./max(U(:,i)));
        end
        
        T = double(Mask);
        for i=1:K
            T(Mask) = U(:,i);
            PC(:,:,i) = T;
        end
        
        pw = reshape(atan2(PC(:,:,1),PC(:,:,2)),dim(1),dim(2));
        Mod =  reshape(sqrt(PC(:,:,1).^2+PC(:,:,2).^2),dim(1),dim(2));
        
    end
    
    varargout{1} = pw;
    varargout{2} = Mod;
    varargout{3} = PC;
    varargout{4} = L;
    varargout{5} = w;
    
catch ME
    throw(ME)
end
