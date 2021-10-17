function [Utrue1,Utrue2] = explicit(Nt,Nx,dx,dt,sigma1,sigma2,r1,r2,A,L,K1)
    % 显格式的求解函数
    
    Utrue1 = zeros(Nt+1,Nx+1);
    Utrue2 = zeros(Nt+1,Nx+1);
    
    alpha = dt/(dx)^2;
    
    beta1 = 1/2-r1/(sigma1^2);
    beta2 = 1/2-r2/(sigma2^2);
    alpha1 = -(r1+A(1,1))-1/(2*sigma1^2)*(r1-sigma1^2/2)^2;
    alpha2 = -(r2+A(2,2))-1/(2*sigma2^2)*(r2-sigma2^2/2)^2;
    xi = alpha2-alpha1;
    eta = beta2-beta1;
    
    for i = 1:Nx+1
        Utrue1(1,i) = exp(-beta1*(-L+(i-1)*dx))*max(K1-exp(-L+(i-1)*dx),0);
        Utrue2(1,i) = exp(-beta2*(-L+(i-1)*dx))*max(K1-exp(-L+(i-1)*dx),0);
    end
    
    for i = 1:Nt+1
        Utrue1(i,Nx+1) = 0;
        Utrue2(i,Nx+1) = 0;
    end
    
    for i = 2:Nt+1
        Utrue1(i,1) = Utrue1(1,1)*exp(-alpha1*(i-1)*dt);
        Utrue2(i,1) = Utrue2(1,1)*exp(-alpha2*(i-1)*dt);
    end
    
    
    for j = 2:Nt+1
        for i = 2:Nx
            Utrue1(j,i) = max(Utrue1(j-1,i)+alpha/2*sigma1^2*(Utrue1(j-1,i+1)-2*Utrue1(j-1,i)+Utrue1(j-1,i-1))-dt*A(1,2)*Utrue2(j-1,i)*exp(xi*(j-2)*dt+eta*(-L+(i-1)*dx)),exp(-(alpha1*(j-1)*dt+beta1*(-L+(i-1)*dx)))*max(K1-exp(-L+(i-1)*dx),0));
            Utrue2(j,i) = max(Utrue2(j-1,i)+alpha/2*sigma2^2*(Utrue2(j-1,i+1)-2*Utrue2(j-1,i)+Utrue2(j-1,i-1))-dt*A(2,1)*Utrue1(j-1,i)*exp(-(xi*(j-2)*dt+eta*(-L+(i-1)*dx))),exp(-(alpha2*(j-1)*dt+beta2*(-L+(i-1)*dx)))*max(K1-exp(-L+(i-1)*dx),0));
        end
    end
    
    for j = 1:Nt+1
        for i = 1:Nx+1
            Utrue1(j,i) = exp(alpha1*(j-1)*dt+beta1*((i-1)*dx-L))*Utrue1(j,i);
            Utrue2(j,i) = exp(alpha2*(j-1)*dt+beta2*((i-1)*dx-L))*Utrue2(j,i);
        end
    end
    
    
end