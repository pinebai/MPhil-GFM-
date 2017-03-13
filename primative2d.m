function [ rho,u,v,p,e ] = primative2d( U )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
gamm = 1.4;
rho = U(:,:,1);
u = U(:,:,2)./U(:,:,1);
v = U(:,:,3)./U(:,:,1);
p = (U(:,:,4) - 0.5*rho.*(u.^2 + v.^2))*(gamm-1);
e = p./((gamm-1).*rho);


end

