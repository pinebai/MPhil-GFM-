clear all
close all
Lx = 1;
Ly = 1;
dx = 0.01;
[xx,yy] = meshgrid([0:dx:Lx],[0:dx:Ly]);
U = zeros([size(xx) 4]);
gamm = 1.4;
t = 0;
CFL = 0.9;
N = size(xx,2);
M = size(xx,1);
w = 0;


rho = 1 - (0.875/2)*(1+erf(1000*(xx-0.5)));
p = 1 - (0.9/2)*(1+erf(1000*(xx-0.5)));
u = xx*0;
v = xx*0;

U(:,:,1) = rho;
U(:,:,2) = u.*rho;
U(:,:,3) = v.*rho;
U(:,:,4) = 0.5*rho.*(u.^2+v.^2) + p/(gamm-1);
t = 0;

while t<0.25
    [U,dt] = SLICstep(U,dx,CFL,w,M,N);
    %[U,dt] = FORCEstep(U,dx,CFL);
    %boundary conditions:
    U(:,end,:) = U(:,end-1,:);
    U(end,:,:) = U(end-1,:,:);
    U(1,:,:) = U(2,:,:);
    U(:,1,:) = U(:,2,:);

    
    t = t+dt
    [ rho,u,v,p,e ] = primative2d( U );
    plot(p(50,:));
    pause(0.001);
    
end





