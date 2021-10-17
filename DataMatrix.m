function [Utrue1,Utrue2] = DataMatrix(Nt,Nx,dx,dt,sigma1,sigma2,r1,r2,gamma,mu,eps,rho0,A,x0,t0,K)
    % 求解出分划网格各点V函数对应的值，存放进Utrue矩阵

    % Nt,Nx分别是时间和空间方向的分段段数
    beta1 = 1/2-r1/(sigma1^2);
    beta2 = 1/2-r2/(sigma2^2);
    
    % DM矩阵用于存放未还原的数值，而Utrue存放还原后的数值
    DM1 = zeros(Nt+1,Nx+1);
    DM2 = zeros(Nt+1,Nx+1);
    Utrue1 = zeros(Nt+1,Nx+1);
    Utrue2 = zeros(Nt+1,Nx+1);
    
    L = XspaceBoundry(sigma1,sigma2,r1,r2,1,K);
    x = linspace(-L,L,Nx+1);
    s = exp(x);
    % 给矩阵左、右、上方三面附上初边值
    % 上方(不包含最右边一个结点)
    DM1(1,1:end-1) = exp(-beta1*x(1:end-1)).*max(K-exp(x(1:end-1)),0);
    DM2(1,1:end-1) = exp(-beta2*x(1:end-1)).*max(K-exp(x(1:end-1)),0);
    % 左方和右方
    for i = 1:Nt+1
        DM1(i,1) = DM1(1,1);
        DM2(i,1) = DM2(1,1);
        DM1(i,end) = 0;
        DM2(i,end) = 0;
    end
    
    % 再给矩阵其余部分赋值
    for i = 2:Nt+1
        W1 = DM1(i-1,2:Nx);
        W2 = DM2(i-1,2:Nx);
        W = [W1 W2].';
        [M,F,G,beta1,beta2,alpha1,alpha2] = AssembleMatrix(i-2,Nx,dx,dt,sigma1,sigma2,r1,r2,A,x0,t0,W,L,K);
        W = W-G;
        Phi_end = PCM(Nx,eps,gamma,mu,rho0,M,F,W);
        DM1(i,2:Nx) = Phi_end(1:Nx-1).'+G(1:Nx-1).';
        DM2(i,2:Nx) = Phi_end(Nx:2*Nx-2).'+G(Nx:2*Nx-2).';
        [V1,V2] = RevValue(i-2,Nx,Phi_end,G,beta1,beta2,alpha1,alpha2,t0,x0,dt,dx);
        Utrue1(i,2:Nx) = V1.';
        Utrue2(i,2:Nx) = V2.';
        
    end
    
    % 还原矩阵边界
    % 还原边界
    for i = 1:Nx
        Utrue1(1,i) = g(s(i),K);
        Utrue2(1,i) = g(s(i),K);
    end
    % 右侧
    for i = 1:Nt+1
        Utrue1(i,Nx+1) = 0;
        Utrue2(i,Nx+1) = 0;
    end
    % 左侧
    for i = 2:Nt+1
        Utrue1(i,1) = Utrue1(1,1);
        Utrue2(i,1) = Utrue2(1,1);
    end
    
end