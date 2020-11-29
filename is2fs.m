function [varargout]=is2fs(varargin);

% IS2FS  Map from image space to feature space.
% CMP Vision Algorithms visionbook@cmp.felk.cvut.cz
% An auxiliary function which takes the data stored as 2D images
% and reshapes it to 1D vector. (Used for PCA demonstration)
% Usage: A = is2fs(X,Xm,U)
% Inputs:
% X    [M x N]  Matrix with the resulting image data. Each column in X
%             is a corresponding vector representing a reconstructed image.
% Xm   [M x 1]  Mean image. 
% U    [M  x K]  Basis vectors of the feature space. 
% Mask [M x N] Processing Mask

% Outputs:
% A  [K  x N]  Coefficients vectors. Each column of 
%              A is a vector in the feature space. 

% Courtesy A. Leonardis, D. Skocaj 
% see http://vicos.fri.uni-lj.si/danijels/downloads


try 
    dim = size(varargin{1});        
    M = dim(1)*dim(2);
    N = dim(3);
    K = varargin{4};    
    
    X  = reshape(varargin{1},M,N);
    Xm = reshape(varargin{2},M,1);
    U  = reshape(varargin{3},M,K);
    
    if Xm==0 Xm=zeros(M,1); end
    A=U'*(X-repmat(Xm,1,N));
    
    varargout{1} = A;
    
catch ME        
    throw(ME)
end
