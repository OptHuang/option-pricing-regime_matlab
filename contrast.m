function contrast
    % 对比不同模型下解的曲线图

    [~,L,Nt,Nx,dx,dt,sigma1,sigma2,r1,r2,gamma,mu,eps,rho0,A,x0,t0,K] = ParaImput();
    
    x = linspace(-L,L,Nx+1);
    s = exp(x);

    Utr1 = zeros(1,Nx+1);
    Utr2 = zeros(1,Nx+1);
    
    [Uxian,~] = explicit(Nt,Nx,dx,dt,sigma1,sigma2,r1,r2,A,L,K);
    
    [Upcm,~] = DataMatrix(Nt,Nx,dx,dt,sigma1,sigma2,r1,r2,gamma,mu,eps,rho0,A,x0,t0,K);
    
    Upx1 = Uxian(end,:);
    Up1 = Upcm(end,:);
    
    for j = 1:Nx+1
        Utr1(j) = TTtransfer(K,Nt,s(j),1,r1,r2,sigma1,sigma2,-A);
        Utr2(j) = TTtransfer(K,Nt,s(j),1,r1,r2,sigma1,sigma2,-A);
    end
    
    plot(s,Utr1(1,:),'r')
    axis([0 25 0 9])
    
    hold on;
    plot(s,Up1(1,:),'y')
    axis([0 25 0 9])
    
    hold on;
    plot(s,Upx1(1,:),'g')
    axis([0 25 0 9])
    
end