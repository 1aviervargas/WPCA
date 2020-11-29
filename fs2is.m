function [varargout]=fs2is(varargin);
% FS2IS  Map from feature space to image space.
% CMP Vision Algorithms http://visionbook.felk.cvut.cz
% An auxiliary function which takes the data stored as a 1D
% vector and reshapes it to a 2D image. (Used for PCA demonstration)
% Usage: X = fs2is(A,Xm,U)
% Inputs:
%   A  [K x N]  Coefficients vectors. Each column of A
%               is a vector in the feature space.
%   Xm  [M x 1]  Mean.
%   U  [M x K]  Basis vectors of the feature space.
% Outputs:
%   X  [M x N]  Matrix with resulting image data. Each column in X
%               is a corresponding vector representing a reconstructed image.

% Courtesy A. Leonardis, D. Skocaj
% see http://vicos.fri.uni-lj.si/danijels/downloads

try
    A = varargin{1};
    Xm = varargin{2};
    U = varargin{3};
    
    [K N] = size(A);        
    nCols = size(Xm,1);
    nRows = size(Xm,2);
    M = nCols*nRows;            
        
    Xm = reshape(Xm,M,1);
    U  = reshape(U,M,K);
    
    N = size( varargin{1}, 2 );
    X = U*A + repmat( Xm, 1, N );
    varargout{1}=reshape(X,nCols,nRows,N);

catch ME        
    throw(ME)
end