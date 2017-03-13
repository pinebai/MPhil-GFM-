function [ U,dt ] = FORCEstep( U,dx,CFL)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
gamm = 1.4;
Uinit = U;

%Using parallel splitting
    [ rho,u,v,p,e ] = primative2d( U );
    %Solving in the x direction

    UL = U(:,1:end-1,:);
    UR = U(:,2:end,:);
    
     a = sqrt(gamm*p./rho);
    cmax = max(max(a + sqrt(u.^2 + v.^2)));
    dt = (CFL*dx)/cmax;
    dt = dt;
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(flux_calc2d(UR,1)-flux_calc2d(UL,1));
    FLW = flux_calc2d(ULW,1);
    
    
    Fforce = 0.25*(flux_calc2d(UL,1) + 2*FLW + flux_calc2d(UR,1) - (dx/dt)*(UR - UL));
    Unew = U;
    
    Unew(:,2:end-1,:) = U(:,2:end-1,:) -(dt/dx)*(Fforce(:,2:end,:)- Fforce(:,1:end-1,:));
    U = Unew;
    
    U(:,end,:) = U(:,end-1,:);
    U(end,:,:) = U(end-1,:,:);
    U(1,:,:) = U(2,:,:);
    U(:,1,:) = U(:,2,:);
    
    %Solving in the y direction
    
    UL = U(1:end-1,:,:);
    UR = U(2:end,:,:);
   
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(flux_calc2d(UR,2)-flux_calc2d(UL,2));
    FLW = flux_calc2d(ULW,2);
    
    
    Fforce = 0.25*(flux_calc2d(UL,2) + 2*FLW + flux_calc2d(UR,2) - (dx/dt)*(UR - UL));
    Unew = U;
    
    Unew(2:end-1,:,:) = U(2:end-1,:,:) -(dt/dx)*(Fforce(2:end,:,:)- Fforce(1:end-1,:,:));
    
    
    U = Unew;
    
    U(:,end,:) = U(:,end-1,:);
    U(end,:,:) = U(end-1,:,:);
    U(1,:,:) = U(2,:,:);
    U(:,1,:) = U(:,2,:);
    
    
     UHS = 0.5*U;
%     
    %secondary splitting:
    
    %Solving in the y direction from initial conditions
    U = Uinit;
   
    UL = U(1:end-1,:,:);
    UR = U(2:end,:,:);
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(flux_calc2d(UR,2)-flux_calc2d(UL,2));
    FLW = flux_calc2d(ULW,2);
    
    
    Fforce = 0.25*(flux_calc2d(UL,2) + 2*FLW + flux_calc2d(UR,2) - (dx/dt)*(UR - UL));
    Unew = U;
    
    Unew(2:end-1,:,:) = U(2:end-1,:,:) -(dt/dx)*(Fforce(2:end,:,:)- Fforce(1:end-1,:,:));
    
    
    U = Unew;
    
    U(:,end,:) = U(:,end-1,:);
    U(end,:,:) = U(end-1,:,:);
    U(1,:,:) = U(2,:,:);
    U(:,1,:) = U(:,2,:);
    
    
    %Solving in x:
    
    UL = U(:,1:end-1,:);
    UR = U(:,2:end,:);
    
    
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(flux_calc2d(UR,1)-flux_calc2d(UL,1));
    FLW = flux_calc2d(ULW,1);
    
    
    Fforce = 0.25*(flux_calc2d(UL,1) + 2*FLW + flux_calc2d(UR,1) - (dx/dt)*(UR - UL));
    Unew = U;
    
    Unew(:,2:end-1,:) = U(:,2:end-1,:) -(dt/dx)*(Fforce(:,2:end,:)- Fforce(:,1:end-1,:));
    U = Unew;
    
    U(:,end,:) = U(:,end-1,:);
    U(end,:,:) = U(end-1,:,:);
    U(1,:,:) = U(2,:,:);
    U(:,1,:) = U(:,2,:);
    
    U = 0.5*U + UHS;

    
end

