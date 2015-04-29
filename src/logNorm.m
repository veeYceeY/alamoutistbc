function [ out ] = logNorm( m,v,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));

out = lognrnd(mu,sigma,1,n);


end

