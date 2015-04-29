function [ gn ] = randGamma( a,b,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for k=1:n
    gn(k)=gamrand(a,b);
end

end

