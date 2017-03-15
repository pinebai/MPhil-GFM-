clear all
close all
Lx = 1;
Ly = 1;
dx = 0.005/sqrt(2);
[xx,yy] = meshgrid([0:dx:Lx],[0:dx:Ly]);
U = zeros([size(xx) 4]);
gamm = 1.4;
t = 0;
CFL = 0.9;
N = size(xx,2);
M = size(xx,1);
w = 0;
Us = [0 0];
b=1;
Test_number = 2;


%       |  rho  |
%U =    | rho u |
%       | rho v |
%       |   E   |


%Sod test x direction

rho = 1 - (0.875/2)*(1+erf(-1000*(1-yy-xx)));
p = 1 - (0.9/2)*(1+erf(-1000*(1-yy-xx)));
u = xx*0;
v = xx*0;

U(:,:,1) = rho;
U(:,:,2) = u.*rho;
U(:,:,3) = v.*rho;
U(:,:,4) = 0.5*rho.*(u.^2+v.^2) + p/(gamm-1);
t = 0;
%phi = sign(abs(real(sqrt(0.05^2 - ((xx-0.3 - Us(1)*t).^2 + (yy-0.5 -Us(2)*t).^2))))-0.001);

phi = 0*xx-1;
%phi(40:60,:) = 1;

phi_s = phi;

for i = 1:10
    phi_s(2:end-1,2:end-1) = phi_s(2:end-1,2:end-1) + 0.1*(phi_s(2:end-1,1:end-2) + phi_s(2:end-1,3:end) + phi_s(1:end-2,2:end-1) + phi_s(3:end,2:end-1) -4*phi(2:end-1,2:end-1));
end
% 
% 
% 
phix = phi*0;
phiy = phi*0;
phix(2:end-1,2:end-1) = (phi_s(2:end-1,3:end)-phi_s(2:end-1,1:end-2))/(2*dx);
phiy(2:end-1,2:end-1) = (phi_s(3:end,2:end-1)-phi_s(1:end-2,2:end-1))/(2*dx);
%Producing a smoothed phi_z for normal calculations:
nx = phix./sqrt(phix.^2 + phiy.^2);
ny = phiy./sqrt(phix.^2 + phiy.^2);
nx(isnan(nx)==1)=0;
ny(isnan(ny)==1)=0;



while t<0.25
    [U,dt] = SLICstep(U,dx,CFL,w,M,N);
    %boundary conditions:
    U(:,end,:) = U(:,end-1,:);
    U(end,:,:) = U(end-1,:,:);
    U(1,:,:) = U(2,:,:);
    U(:,1,:) = U(:,2,:);

    
    t = t+dt
    %axis([0 1 0 1 0 1]) 
    [ rho,u,v,p,e ] = primative2d( U );
    
    %creating extrapolated values:
    rho_ext = rho;
    u_ext = u;
    v_ext = v;
    p_ext = p;
    dtau = dt/10;
    
%     phi = sign(abs(real(sqrt(0.05^2 - ((xx-0.4 - Us(1)*t).^2 + (yy-0.5 -Us(2)*t).^2))))-0.001);
%     
%     phi_s = phi;
%     
%     for i = 1:10
%         phi_s(2:end-1,2:end-1) = phi_s(2:end-1,2:end-1) + 0.1*(phi_s(2:end-1,1:end-2) + phi_s(2:end-1,3:end) + phi_s(1:end-2,2:end-1) + phi_s(3:end,2:end-1) -4*phi(2:end-1,2:end-1));
%     end
%     
%     
%     
%     phix = phi*0;
%     phiy = phi*0;
%     phix(2:end-1,2:end-1) = (phi_s(2:end-1,3:end)-phi_s(2:end-1,1:end-2))/(2*dx);
%     phiy(2:end-1,2:end-1) = (phi_s(3:end,2:end-1)-phi_s(1:end-2,2:end-1))/(2*dx);
%     %Producing a smoothed phi_z for normal calculations:
%     nx = phix./sqrt(phix.^2 + phiy.^2);
%     ny = phiy./sqrt(phix.^2 + phiy.^2);
%     nx(isnan(nx)==1)=0;
%     ny(isnan(ny)==1)=0;
    
    
    for tau = 1:100
        for i = 2:N-1
            for j = 2:M-1
                if phi_s(j,i)>=0 %ie we are in the ghost fluid region
                    disp('error')
                    if nx(j,i)>0 %see appendix A.4 in Fedkiw et al.
                        %calculating upwinded derivatives
                        drhodx = (rho(j,i) - rho(j,i-1))/dx;
                        dpdx = (p(j,i) - p(j,i-1))/dx;
                        %If the velocity derivative uses a value from
                        %outside of the ghost cell region, the tangential
                        %component of velocity is interpolated normally,
                        %but the normal component of velocity must be
                        %reflected before being interpolated
                        
                        if phi_s(j,i-1)>=0 %ie both cells are ghost cells
                            dudx = (u(j,i) - u(j,i-1))/dx;
                            dvdx = (v(j,i) - v(j,i-1))/dx;
                        else
                            %Initially finding the reflected tangential component of
                            %velocity of u(j,i-1):
                            urel = u(j,i-1) - Us(1);
                            vrel = v(j,i-1) - Us(2);
                            uref = b*(urel - 2*(urel*nx(j,i) + vrel*ny(j,i))*nx(j,i));
                            vref = b*(vrel - 2*(urel*nx(j,i) + vrel*ny(j,i))*ny(j,i));
                            dudx = (u(j,i) - uref)/dx;
                            dvdx = (v(j,i) - vref)/dx;
                        end
                        
                        
                    elseif nx(j,i)<0
                        drhodx = (rho(j,i+1) - rho(j,i))/dx;
                        dpdx = (p(j,i+1) - p(j,i))/dx;
                        if phi_s(j,i+1)>=0 %ie both cells are ghost cells
                            dudx = (u(j,i+1) - u(j,i))/dx;
                            dvdx = (v(j,i+1) - v(j,i))/dx;
                        else
                            %Initially finding the reflected tangential component of
                            %velocity of u(j,i-1):
                            urel = u(j,i+1) - Us(1);
                            vrel = v(j,i+1) - Us(2);
                            uref = b*(urel - 2*(urel*nx(j,i) + vrel*ny(j,i))*nx(j,i));
                            vref = b*(vrel - 2*(urel*nx(j,i) + vrel*ny(j,i))*ny(j,i));
                            dudx = (uref - u(j,i))/dx;
                            dvdx = (vref - v(j,i))/dx;
                        end
                    else
                        dpdx = 0;
                        dudx = 0;
                        drhodx = 0;
                        dvdx = 0;
                    end
                    
                    if ny(j,i)>0
                        drhody = (rho(j,i) -rho(j-1,i))/dx;
                        dpdy = (p(j,i)-p(j-1,i))/dx;
                        if phi_s(j-1,i)>=0
                            dudy = (u(j,i)-u(j-1,i))/dx;
                            dvdy = (v(j,i)-v(j-1,i))/dx;
                        else
                            urel = u(j-1,i) - Us(1);
                            vrel = v(j-1,i) - Us(2);
                            uref = b*(urel - 2*(urel*nx(j,i) + vrel*ny(j,i))*nx(j,i));
                            vref = b*(vrel - 2*(urel*nx(j,i) + vrel*ny(j,i))*ny(j,i));
                            dudy = (u(j,i)-uref)/dx;
                            dvdy = (v(j,i)-vref)/dx;

                        end
                    elseif ny(j,i)<0
                        drhody = (rho(j+1,i)-rho(j,i))/dx;
                        dpdy = (p(j+1,i)-p(j,i))/dx;
                        if phi_s(j+1,i)>=0
                            dudy = (u(j+1,i)-u(j,i))/dx;
                            dvdy = (v(j+1,i)-v(j,i))/dx;
                        else
                            urel = u(j+1,i) - Us(1);
                            vrel = v(j+1,i) - Us(2);
                            uref = b*(urel - 2*(urel*nx(j,i) + vrel*ny(j,i))*nx(j,i));
                            vref = b*(vrel - 2*(urel*nx(j,i) + vrel*ny(j,i))*ny(j,i));
                            dudy = (uref - u(j,i))/dx;
                            dvdy = (vref - v(j,i))/dx;
                        end
                    else 
                        dudy = 0;
                        dvdy = 0;
                        dpdy = 0;
                        drhody = 0;
                    end
                    
                    u_ext(j,i) = u(j,i) - dtau*(dudy.*ny(j,i) + dudx.*nx(j,i));
                    v_ext(j,i) = v(j,i) - dtau*(dvdy.*ny(j,i) + dvdx.*nx(j,i));
                    p_ext(j,i) = p(j,i) - dtau*(dpdy.*ny(j,i) + dpdx.*nx(j,i));
                    rho_ext(j,i) = rho(j,i) - dtau*(drhody*ny(j,i) + drhodx*nx(j,i)); 
                end
            end
        end

        u = u_ext;
        v = v_ext;
        p = p_ext;
%         rho = rho_ext;
%         clf
%         hold on
%         surf(p)
%         surf(phi)
%         axis([0 100 0 100 0 1])
%         view(-37.5,60);
%         pause(0.1)
        
        
    end
    U(:,:,1) = rho;
    U(:,:,2) = u.*rho;
    U(:,:,3) = v.*rho;
    U(:,:,4) = 0.5*rho.*(u.^2+v.^2) + p/(gamm-1);
    
                       
                    
                    
                    
                    
                   
%         u = u_ext;
%         v = v_ext;
%         p = p_ext;
%         rho = rho_ext;
%         clf
%         hold on
%         %surf(p)
%         %plot(p(50,:))
%         %surf(phi*2.6)
%        axis([0 100 0 100 0 2.5])
%         view(-37.5,60);
%         pause(0.1)
%surf(p.*(sign(-phi)+1)/2)

%imagesc(p.*(sign(-phi)+1)/2);
%mock schlieren:

grad_rho = abs((rho(3:end,2:end-1) - rho(1:end-2,2:end-1)))/(2*dx) + abs(rho(2:end-1,3:end)-rho(2:end-1,1:end-2))/(2*dx);
schil = exp(-20*abs(grad_rho)./(1000*rho(2:end-1,2:end-1)));
schil = schil.*(sign(-phi(2:end-1,2:end-1))+1)/2;
schil(schil==0)=1;
%surf(schil);
%image(50*schil);
%imagesc(schil.*(sign(-phi(2:end-1,2:end-1))+1)/2);
%imagesc(schil);
% plot(rho(:,37));
% hold on
% plot(rho(:,65));
% %plot(rho(:,32),'x');
% hold off
%surf(p);
%quiver(u,v);
%view(60,60);
schil_pad = ones(M,N);
%imagesc(-log(schil));
%plot(p(50,:));
schil_pad(2:end-1,2:end-1) = schil;
pause(0.0001);
end

%centerline plotter:

for i = 1:min([N,M])
    pdiag(i) = p(i,i);
    xdiag(i) = sqrt(2)*xx(1,i)-(sqrt(2)-1)/2;
end
        


