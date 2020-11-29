function varargout = adjustPiston(varargin)

%ADJUSTPISTON Adjust the piston term between two wrapped phases
%and obtains also the rms error between them
%
% Usage: [pwAdjusted rmsAdjusted] = adjustPiston(pwRef,pw)
%
% Inputs:
%
%   pwRef [NRows x NCols] Reference wrapped phase
%
%   pw    [NRows x NCols] Wrapped phase to adjust to pwRef
%       and columns of the fringe pattern   
%
%   Mask [NRows x NCols] Processing mask
%
% Outputs:
%   pwAdjust  [NRows x NCols]  Adjusted wrapped phase.
%
%   rmsAdjust [1 x 1 value]  Obtained rms error between pwAdjust and pwRef
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

try
    
    pwRef = varargin{1};
    pw = varargin{2};
    
    if (length(varargin) < 3 )
        Mask = ones(size(pw))>0.5;
    else
        Mask = varargin{3};
    end
    
    pwRef = mod(pwRef,2*pi);
    pw = mod(pw,2*pi);
    
    index =1;
    range = 0:0.1:2*pi;
    rmsVectorPlus=zeros(1,length(range));
    rmsVectorMinus=zeros(1,length(range));
    
    for i=range        
        %rmsVectorPlus(index)=sqrt(mean2(abs(pwRef.*Mask-mod(pw+i,2*pi).*Mask)));
        %rmsVectorMinus(index)=sqrt(mean2(abs(pwRef.*Mask-mod(-pw+i,2*pi).*Mask)));
        rmsVectorPlus(index)=sqrt(mean2(abs(pwRef(Mask)-mod(pw(Mask)+i,2*pi))));
        rmsVectorMinus(index)=sqrt(mean2(abs(pwRef(Mask)-mod(-pw(Mask)+i,2*pi))));
        index=index+1;
    end
    
    [rmsPlus, indxPlus] = min(rmsVectorPlus);
    [rmsMinus, indxMinus] = min(rmsVectorMinus);
    
    sign = nan;
    rmsAdjusted = nan;
    
    if (rmsPlus < rmsMinus )
        ps = range(indxPlus);
        rmsAdjusted = rmsPlus;
        sign = 1;
    else
        ps = range(indxMinus);
        sign = -1;
        rmsAdjusted = rmsMinus;
    end
    
    pwAdjusted = mod(sign*pw+ps,2*pi);
    
    varargout{1} = pwAdjusted;
    varargout{2} = rmsAdjusted;
    
catch ME
    throw(ME)
end
    