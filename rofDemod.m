function varargout = rofDemod(varargin)

%ROFDEMOD Fringe pattern demodulation method from two phase-shifted
%interferograms. This function implements the method presented in [1]. 
%The method is based on first obtaining the direction map by a regularized
%optical flow algorithm [2]. Then the Spiral Phase Transform [3] is applied
%to one of the fringe patterns and the wrapped phase it is obtained.
%
% Usage: [pw,Mod] = rofDemod(I1,I2,K,epsilon,lambda,Mask)
%
% Inputs:
%
%   I1  [NRows x NCols] Image map with NRows and NCols the number of rows 
%       and columns of the fringe pattern
%   I2  [NRows x NCols] Image map with NRows and NCols the number of rows 
%       and columns of the fringe pattern   
%   K   [1x1 double value] Max number of iterations of the regularization
%       process. (Default ~ 5-10)
%   elpsilon [1x1 double value] Stopping criteria of the regularization
%       process. (Default ~ 0.3)
%   lambda [1x1 double value] Regularization parameter (Default ~ 10)
%
%   Mask [NRows x NCols] Processing mask
%
% Outputs:
%   pw  [NRows x NCols]  Wrapped modulating phase.
%   Mod [NRows x NCols]  Modulation term
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
%
%[1] J. Vargas, J. A Quiroga, C.O.S. Sorzano, J. C. Estrada, J. M. Carazo,
%    "Two-step interferometry by a regularized optical flow algorithm",
%    sent for publication to Optics Letters (2011)
%
%[2] B. K. Horn, "Determining Optical Flow", "Artificial Intelligence", 81,
%   (1981). (http://www.sciencedirect.com/science/article/pii/0004370281900242)
%
%[3] Kieran G. Larkin, "Natural demodulation of two-dimensional fringe 
%    patterns. II. Stationary phase analysis of the spiral phase quadrature 
%    transform", JOSA A, Vol. 18, Issue 8, pp. 1871-1881 (2001)

try
    
    switch length(varargin)
        case 1
            error('At least are necessary to images to demodulated them by the regularized optical flow method');
        case 2
            K=10;
            epsilon=0.3;
            lambda=15;
            Mask = ones(size(varargin{1})) > 0.5;
        case 3
            K = varargin{3};
            epsilon=0.3;
            lambda=15;
            Mask = ones(size(varargin{1})) > 0.5;
        case 4
            K = varargin{3};
            epsilon=varargin{4};
            lambda=15;
            Mask = ones(size(varargin{1})) > 0.5;
        case 5
            K = varargin{3};
            epsilon=varargin{4};
            lambda=varargin{5};
            Mask = ones(size(varargin{1})) > 0.5;
        case 6
            K = varargin{3};
            epsilon=varargin{4};
            lambda=varargin{5};
            Mask = varargin{6};
    end
    
    I1 = varargin{1};
    I2 = varargin{2};
    
    [U,V] = OptFlowGS(I1,I2,K,epsilon,lambda);
    D = medfilt2(atan2(-V,U));
    
    % DC filtering
    [NR, NC]=size(I1);
    [u,v]=meshgrid(1:NC, 1:NR);
    u0=floor(NC/2)+1; v0=floor(NR/2)+1;
    u=u-u0; v=v-v0;
       
    H=1-exp(-(u.^2+v.^2)/(2*0.05^2)); %Gaussian DC filter with sigma=0.05
    C=fft2(I1);
    CH=C.*ifftshift(H);
    ch=ifft2(CH);
    
    %compute quadrature using the Spiral PHase Transform
    sd=(SPHT(real(ch)));
    s=-1i.*exp(-1i*D).*sd;
    
    pROF = atan2(real(s), real(ch));
    pROF = pROF.*Mask;
    
    Mod = sqrt(real(s).^2+real(ch).^2);
    Mod = (Mod / max(Mod(:))).*Mask;
    
    varargout{1} = pROF;
    varargout{2} = Mod;
    
catch ME   
    throw(ME)
end    