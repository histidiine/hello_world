function [output] = Get_MinSize(varargin)

T_size= zeros(nargin,1);
for ii = 1:nargin;
    T_size(ii,1)=length(varargin{ii});
    ii = ii+1;
end

output = min(T_size);