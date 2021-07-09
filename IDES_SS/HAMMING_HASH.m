function [out] = HAMMING_HASH(perm_key,w)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

out=[];
for i=1:length(perm_key)
    xused=perm_key(i);
    phi=w(xused);
    out=[out phi];
    
end


end

