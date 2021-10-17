function Err
    % 求解误差和收敛阶
    
    [T,L,Nt,Nx,dx,dt,sigma1,sigma2,r1,r2,gamma,mu,eps,rho0,A,x0,t0,K] = ParaImput();

    Norm = zeros(1,5);
    
    for i = 1:5
    x = linspace(-L,L,i*Nx+1);
    s = exp(x);

    Utr1 = zeros(1,i*Nx+1);
    Utr2 = zeros(1,i*Nx+1);
    [Upcm1,Upcm2] = DataMatrix(Nt,i*Nx,dx/i,dt,sigma1,sigma2,r1,r2,gamma,mu,eps,rho0,A,x0,t0,K);
    
    Up1 = Upcm1(Nt+1,:);
    Up2 = Upcm2(Nt+1,:);

    for j = 1:i*Nx+1
        Utr1(j) = TTtransfer(K,Nt,s(j),T,r1,r2,sigma1,sigma2,-A);
        Utr2(j) = TTtransfer(K,Nt,s(j),T,r1,r2,sigma1,sigma2,-A);
    end    
    intF = 0.5*(Utr1(1:end-1)+Utr1(2:end)-Up1(1:end-1)-Up1(2:end)); 
    Norm(i) = sqrt(intF*intF'*2*L/(Nx*i));
    end
    
    N0 = [50 100 150 200 250];
    loglog(2*L./N0,Norm,'-+')
   
%     N = [60 70 80 90 100];
%     Err1
%     
%     loglog((2*L)./N,Err1,'-+')
%     hold on
%     loglog((2*L)./N,Err2,'-+')
%         
end