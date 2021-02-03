function varargout = gsDemod(varargin)
    
%GSDEMOD Fringe pattern demodulation method from two phase-shifted
%interferograms. This function implements the method presented in [1]. 
%The method is based on orthonormalizing the fringe patterns using the 
%Gram-Schmidt orthogonalization method. If there is more that one fringe,
%this process is equivalent to demodulate the fringe patterns.
%
% Usage: [pw,Mod] = rofDemod(I1,I2,Mask)
%
% Inputs:
%
%   I1  [NRows x NCols] Image map with NRows and NCols the number of rows 
%       and columns of the fringe pattern
%   I2  [NRows x NCols] Image map with NRows and NCols the number of rows 
%       and columns of the fringe pattern   
%   Mask [NRows x NCols] Processing mask
%
% Outputs:
%   pw  [NRows x NCols]  Wrapped modulating phase.
%   Mod [NRows x NCols]  Modulation term
%
%   Javier Vargas, 
%   24/08/2011 
%   copyright @2011 
%   Biocomputing Unit, Centro Nacional de Biotecnología (CSIC)
%   http://biocomp.cnb.csic.es/
%   $ Revision: 1.0.0.0 $
%   $ Date: 25/10/10 $
%   $ Revision: 1.0.0.0 $
%   $ Date: 24/08/2011 $
%
% THIS CODE IS GIVEN WITHOUT ANY GUARANTY AND WITHOUT ANY SUPPORT
%
% REFERENCES 
%
%[1] J. Vargas, J. A Quiroga, C.O.S. Sorzano, J. C. Estrada, J. M. Carazo,
%    "Two-step demodulation based on Gram-Schmidt orthonormalization method",
%    sent for publication to Optics Letters (2011)

try
    
    
    switch length(varargin)
       
        case 1
            error('At least are necessary to images to demodulated them by the Gram-Schmidt method');
        case 2
            I1 = varargin{1};
            I2 = varargin{2};
            Mask = ones(size(varargin{1})) > 0.5;
        case 3
            I1 = varargin{1};
            I2 = varargin{2};
            Mask = varargin{3};
        otherwise
            error('This function requires as arguments two images and a processing mask');
                
    end
    
    norm1 = sqrt(sum(I1(Mask).*I1(Mask)));
    I1 = I1./norm1;
    
    proj = sum((I1(Mask).*I1(Mask)))./sum(I1(Mask).*I2(Mask));
    
    I2 = I1-proj*I2;
    norm2 = sqrt(sum(I2(Mask).*I2(Mask)));
    I2 = I2./norm2;
    
    pw = atan2(-I2,I1);

    M = abs(I2)+abs(I1);
    PC(:,:,1) = -I2;
    PC(:,:,2) = I1;
    
    varargout{1}=pw;
    varargout{2}=M;
    varargout{3}=PC;
    
catch ME
    throw(ME)
end
