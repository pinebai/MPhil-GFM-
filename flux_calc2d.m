function [ F ] = flux_calc2d( U,dim )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
gamm = 1.4;
rho = U(:,:,1);
u = U(:,:,2)./U(:,:,1);
v = U(:,:,3)./U(:,:,1);
p = (U(:,:,4) - 0.5*rho.*(u.^2+v.^2))*(gamm-1);
if dim == 1
    F(:,:,1) = rho.*u;
    F(:,:,2) = rho.*u.^2 +p;
    F(:,:,3) = rho.*u.*v;
    F(:,:,4) = u.*(U(:,:,4)+p);
    %F(:,:,5) = u.*rho.*U(:,:,5);
else
    F(:,:,1) = rho.*v;
    F(:,:,3) = rho.*v.^2 +p;
    F(:,:,2) = rho.*u.*v;
    F(:,:,4) = v.*(U(:,:,4)+p);
    %F(:,:,5) = v.*rho.*U(:,:,5);
end

end

