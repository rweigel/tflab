function PE = pe_nonflag(y,f,varargin)
%PE_NONFLAG Prediction efficiency between non-FLAG elements
% 
%   PE_NONFLAG(A,P) returns 1-ARV_NONFLAG(A,P)
%
%   PE_NONFLAG(A,P,FLAG,EXACT) See IS_FLAG for a definition of FLAG and COND.
%  
%   If all elements along dimension DIM are equal to FLAG, a fill of FLAG is
%   used.
%   
%   See also ARV, MSE, *_NONFLAG.

% R.S. Weigel, 04/02/2004.

PE = 1-arv_nonflag(y,f,varargin{:});  