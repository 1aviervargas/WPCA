function varargout = kreisDemod(varargin)

%KREISDEMOD Fringe pattern demodulation method from two phase-shifted
%interferograms using the Kreis method [1]
%
% Usage: [pw,Mod] = kreisDemod(I1,I2,Mask)
%
% Inputs:
%
%   I1  [NRows x NCols] Image map with NRows and NCols the number of rows
%       and columns of the fringe pattern
%   Mask [NRows x NCols] Processing mask
%
% Outputs:
%   pw  [NRows x NCols]  Wrapped modulating phase.
%
%   Javier Vargas,
%   06/06/2011
%   copyright @2011
%   Biocomputing Unit, Centro Nacional de Biotecnolog�a (CSIC)
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
%[1] T. M. Kreis, P. P Jueptner, "Fourier transform evaluation of 
%    interference patterns: demodulation and sign ambiguity", Proc of SPIE,
%    vol. 1553 (1992)

try
    
    I1 = varargin{1};
    I2 = varargin{2};
    
    if (length(varargin) < 3 )
        Mask = ones(size(I1)) > 0.5;
    else
        Mask = varargin{3};
    end
    
    %Compute the orientation terms
    [fx1 fy1] = gradient(I1);
    [fx2 fy2] = gradient(I2);
    o1 = atan2(-fy1,fx1);
    
    [NR, NC]=size(I1);
    [u,v]=meshgrid(1:NC, 1:NR);
    u0=floor(NC/2)+1; v0=floor(NR/2)+1;
    u=u-u0; v=v-v0;
    
    %Gaussian DC filter with sigma=0.1 to substract the DC component
    H=1-exp(-(u.^2+v.^2)/(2*0.1^2));
    C=fft2(I1);
    CH=C.*ifftshift(H);
    ch1=ifft2(CH);
    
    C=fft2(I2);
    CH=C.*ifftshift(H);
    ch2=ifft2(CH);
    
    %compute quadrature using the Spiral PHase Transform
    s1=(SPHT(real(ch1)));
    s2=(SPHT(real(ch2)));
    
    p1 = atan2(real(-1i.*exp(-1i*o1).*s1), real(ch1));
    p2 = atan2(real(-1i.*exp(-1i*o1).*s2), real(ch2));
    
    %Alpha=mod(p2-p1-pi/2, pi)-pi/2;
    c1 = exp(i*p1);
    c2 = exp(i*p2);
    Alpha = atan2((real(c1).*imag(c2)-imag(c1).*real(c2)),(real(c1).*real(c2)+imag(c1).*imag(c2)) );
    %quadrature sign
    QS=sign(Alpha);
    
    %retrieved phase by the Kresis method
    pK=(QS.*p1).*Mask;
    PC(:,:,1) = -sin(pK);
    PC(:,:,2) = cos(pK);
    Mod = sqrt( (PC(:,:,1)).^2+(PC(:,:,2)).^2);
    
    varargout{1} = pK;
    varargout{2} = Mod;
    varargout{3} = PC;
    
catch ME
    throw(ME)
end
