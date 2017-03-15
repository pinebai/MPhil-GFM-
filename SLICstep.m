function [ U,dt ] = SLICstep( U,dx,CFL,w,M,N)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
gamm = 1.4;
Uinit = U;

%Using parallel splitting, ie: e^(A + B) \approx 0.5*(e^(A)e^(B) +
%e^(B)e^(A).

    [ rho,u,v,p,e ] = primative2d( U );
    a = sqrt(gamm*p./rho);
    cmax = max(max(a + sqrt(u.^2 + v.^2)));
    dt = (CFL*dx)/cmax;
    Beta = 1;
    
    %Solving in the x direction initially
    
    %creating a TVD flux limiter (see Toro for details)
    delta_i2 = U(:,2:end,:)-U(:,1:end-1,:);
    r = delta_i2(:,1:end-1,:)./delta_i2(:,2:end,:);
    r(isnan(r)==1)=0; %removing nan values from 0/0;
    zetasb = zeros(M,N-2); %Creating superbee flux vector
    Beta = 1; %Limiter parameter
    ZetaR = 2*Beta./(1-w +(1+w).*r); %Rightward upwinded flux limiter
    
    %looping over all values of r to create the flux limitng variable:
    k = 1; 
    %k determines which of the primitive variables will be 
    %used for flux limiting, default is k = 1, ie oscillations in the
    %density will be used for the flux limiter calculations
    
    for i = 1:N-2
        for j = 2:M-1
            if r(j,i,k)<=0
                zetasb(j,i) = 0;
            elseif r(j,i,k) <=0.5
                zetasb(j,i) = 2*r(j,i,k);
            elseif r(j,i,k) <= 1;
                zetasb(j,i) = 1;
            else
                zetasb(j,i) = min([r(j,i,k) ZetaR(j,i,k) 2]);
            end
        end
    end
    
   %Applying flux limiter to interpolated values:
   for i = 1:4
       delta_i(:,1:N-2,i) = zetasb(:,1:N-2).*(0.5*(1+w).*delta_i2(:,1:N-2,i) + 0.5*(1-w).*delta_i2(:,2:N-1,i));
   end
    
   %Stepping through half a time step using HLLC fluxes
    ULs = U(:,2:end-1,:) - 0.5*delta_i;
    URs = U(:,2:end-1,:) + 0.5*delta_i;
    
    ULbar = ULs + 0.5*(dt/dx)*(flux_calc2d(ULs,1) - flux_calc2d(URs,1));
    URbar = URs + 0.5*(dt/dx)*(flux_calc2d(ULs,1) - flux_calc2d(URs,1));
    
    UL = URbar(:,1:end-1,:);
    UR = ULbar(:,2:end,:);
    
    %Using force to step through the rest of the time step.
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(flux_calc2d(UR,1)-flux_calc2d(UL,1));
    FLW = flux_calc2d(ULW,1);
    
    
    Fforce = 0.25*(flux_calc2d(UL,1) + 2*FLW + flux_calc2d(UR,1) - (dx/dt)*(UR - UL));
    Unew = U;
    
    Unew(:,3:end-2,:) = U(:,3:end-2,:) -(dt/dx)*(Fforce(:,2:end,:)- Fforce(:,1:end-1,:));
    U = Unew;
    
    U(:,end,:) = U(:,end-2,:);
    U(:,end-1,:) = U(:,end-2,:);
    U(end,:,:) = U(end-2,:,:);
    U(end-1,:,:) = U(end-2,:,:);
    U(1,:,:) = U(3,:,:);
    U(2,:,:) = U(3,:,:);
    U(:,1,:) = U(:,3,:);
    U(:,2,:) = U(:,3,:);
    
    %Solving in the y direction

    delta_i2 = U(2:end,:,:) - U(1:end-1,:,:);
    r = delta_i2(1:end-1,:,:)./delta_i2(2:end,:,:);
    r(isnan(r)==1)=0;
    zetasb = zeros(M-2,N);
    zetaR = 2*Beta./(1-w + (1+w).*r);
    delta_i = [];
    
    %again using density for TVD limiter:
    k = 1;
    for i = 2:N-1
        for j = 1:M-2
            if r(j,i,k)<=0
                zetasb(j,i) = 0;
            elseif r(j,i,k)<=0.5
                zetasb(j,i) = 2*r(j,i,k);
            elseif r(j,i,k)<=1
                zetasb(j,i) = 1;
            else
                zetasb(j,i) = min([r(j,i,k) zetaR(j,i,k) 2]);
            end
        end
    end

    for i = 1:4
        delta_i(1:M-2,:,i) = zetasb(1:M-2,:).*(0.5*(1+w).*delta_i2(1:M-2,:,i) + 0.5*(1-w).*delta_i2(2:M-1,:,i));
    end
    
    ULs = U(2:end-1,:,:) - 0.5*delta_i;
    URs = U(2:end-1,:,:) + 0.5*delta_i;
    
    ULbar = ULs +0.5*(dt/dx)*(flux_calc2d(ULs,2) - flux_calc2d(URs,2));
    URbar = URs + 0.5*(dt/dx)*(flux_calc2d(ULs,2) - flux_calc2d(URs,2));
    
    UL = URbar(1:end-1,:,:);
    UR = ULbar(2:end,:,:);
    
    
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(flux_calc2d(UR,2)-flux_calc2d(UL,2));
    FLW = flux_calc2d(ULW,2);
    
    
    Fforce = 0.25*(flux_calc2d(UL,2) + 2*FLW + flux_calc2d(UR,2) - (dx/dt)*(UR - UL));
    Unew = U;
    
    Unew(3:end-2,:,:) = U(3:end-2,:,:) -(dt/dx)*(Fforce(2:end,:,:)- Fforce(1:end-1,:,:));
    
    
    U = Unew;
    
    U(:,end,:) = U(:,end-2,:);
    U(:,end-1,:) = U(:,end-2,:);
    U(end,:,:) = U(end-2,:,:);
    U(end-1,:,:) = U(end-2,:,:);
    U(1,:,:) = U(3,:,:);
    U(2,:,:) = U(3,:,:);
    U(:,1,:) = U(:,3,:);
    U(:,2,:) = U(:,3,:);
    
    
     UHS = 0.5*U;
     
     U = Uinit;
   delta_i2 = U(2:end,:,:) - U(1:end-1,:,:);
    r = delta_i2(1:end-1,:,:)./delta_i2(2:end,:,:);
    r(isnan(r)==1)=0;
    zetasb = zeros(M-2,N);
    zetaR = 2*Beta./(1-w + (1+w).*r);
    delta_i = [];
    %delta_i = zeros(
    
    %again using density for TVD limiter:
    k = 1;
    for i = 2:N-1
        for j = 1:M-2
            if r(j,i,k)<=0
                zetasb(j,i) = 0;
            elseif r(j,i,k)<=0.5
                zetasb(j,i) = 2*r(j,i,k);
            elseif r(j,i,k)<=1
                zetasb(j,i) = 1;
            else
                zetasb(j,i,k) = min([r(j,i,k) zetaR(j,i,k) 2]);
            end
        end
    end

    for i = 1:4
        delta_i(:,:,i) = zetasb(1:M-2,:).*(0.5*(1+w).*delta_i2(1:M-2,:,i) + 0.5*(1-w).*delta_i2(2:M-1,:,i));
    end
    
    ULs = U(2:end-1,:,:) - 0.5*delta_i;
    URs = U(2:end-1,:,:) + 0.5*delta_i;
    
    ULbar = ULs +0.5*(dt/dx)*(flux_calc2d(ULs,2) - flux_calc2d(URs,2));
    URbar = URs + 0.5*(dt/dx)*(flux_calc2d(ULs,2) - flux_calc2d(URs,2));
    
    UL = URbar(1:end-1,:,:);
    UR = ULbar(2:end,:,:);
    
    
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(flux_calc2d(UR,2)-flux_calc2d(UL,2));
    FLW = flux_calc2d(ULW,2);
    
    
    Fforce = 0.25*(flux_calc2d(UL,2) + 2*FLW + flux_calc2d(UR,2) - (dx/dt)*(UR - UL));
    Unew = U;
    
    Unew(3:end-2,:,:) = U(3:end-2,:,:) -(dt/dx)*(Fforce(2:end,:,:)- Fforce(1:end-1,:,:));
    
    
    U = Unew;
    
    U(:,end,:) = U(:,end-2,:);
    U(:,end-1,:) = U(:,end-2,:);
    U(end,:,:) = U(end-2,:,:);
    U(end-1,:,:) = U(end-2,:,:);
    U(1,:,:) = U(3,:,:);
    U(2,:,:) = U(3,:,:);
    U(:,1,:) = U(:,3,:);
    U(:,2,:) = U(:,3,:);
    
    %Solving in the x direction:
    
    delta_i2 = U(:,2:end,:)-U(:,1:end-1,:);
    r = delta_i2(:,1:end-1,:)./delta_i2(:,2:end,:);
    r(isnan(r)==1)=0; %removing nan values from 0/0;
    zetasb = zeros(M,N-2); %Creating superbee flux vector
    Beta = 1; %Limiter parameter
    ZetaR = 2*Beta./(1-w +(1+w).*r); %Rightward upwinded flux limiter
    delta_i = [];
    %looping over all values of r to create the flux limitng variable:
    k = 1; 
    %k determines which of the primitive variables will be 
    %used for flux limiting, default is k = 1, ie oscillations in the
    %density will be used for the flux limiter calculations
    
    for i = 1:N-2
        for j = 2:M-1
            if r(j,i,k)<=0
                zetasb(j,i) = 0;
            elseif r(j,i,k) <=0.5
                zetasb(j,i) = 2*r(j,i,k);
            elseif r(j,i,k) <= 1;
                zetasb(j,i) = 1;
            else
                zetasb(j,i) = min([r(j,i,k) ZetaR(j,i,k) 2]);
            end
        end
    end
    
   %Applying flux limiter to interpolated values:
   for i = 1:4
       delta_i(:,1:N-2,i) = zetasb(:,1:N-2).*(0.5*(1+w).*delta_i2(:,1:N-2,i) + 0.5*(1-w).*delta_i2(:,2:N-1,i));
   end
    
   %Stepping through half a time step using HLLC fluxes
    ULs = U(:,2:end-1,:) - 0.5*delta_i;
    URs = U(:,2:end-1,:) + 0.5*delta_i;
    
    ULbar = ULs + 0.5*(dt/dx)*(flux_calc2d(ULs,1) - flux_calc2d(URs,1));
    URbar = URs + 0.5*(dt/dx)*(flux_calc2d(ULs,1) - flux_calc2d(URs,1));
    
    UL = URbar(:,1:end-1,:);
    UR = ULbar(:,2:end,:);
    
    %Using force to step through the rest of the time step.
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(flux_calc2d(UR,1)-flux_calc2d(UL,1));
    FLW = flux_calc2d(ULW,1);
    
    
    Fforce = 0.25*(flux_calc2d(UL,1) + 2*FLW + flux_calc2d(UR,1) - (dx/dt)*(UR - UL));
    Unew = U;
    
    Unew(:,3:end-2,:) = U(:,3:end-2,:) -(dt/dx)*(Fforce(:,2:end,:)- Fforce(:,1:end-1,:));
    U = Unew;
    
    U(:,end,:) = U(:,end-2,:);
    U(:,end-1,:) = U(:,end-2,:);
    U(end,:,:) = U(end-2,:,:);
    U(end-1,:,:) = U(end-2,:,:);
    U(1,:,:) = U(3,:,:);
    U(2,:,:) = U(3,:,:);
    U(:,1,:) = U(:,3,:);
    U(:,2,:) = U(:,3,:);
    
    %Combining split equations:
    
    U = 0.5*U + UHS;
 

    
end

