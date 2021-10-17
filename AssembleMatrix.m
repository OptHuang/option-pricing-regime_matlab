function [M,F,G,beta1,beta2,alpha1,alpha2] = AssembleMatrix(n,Nx,dx,dt,sigma1,sigma2,r1,r2,A,x0,t0,W,L,K)
    % 形成矩阵M和F

    % n是现阶段迭代步骤所在的时间方向的层数(n层到n+1层)
    % Nx是空间方向剖分段数，dx是每小段长度，x0是最左侧点
    % dt是时间分划每小段长度，t0是最下方起始点
    % sigma1,sigma2,r1,r2是给定系数,A是给定的转换矩阵
    % W是一个(2Nx-2)*1的列向量，是时间层第n层的数据
    % L是截取边界，K是常参数
    
    % 坐标变换导致的参数变化
    beta1 = 1/2-r1/(sigma1^2);
    beta2 = 1/2-r2/(sigma2^2);
    alpha1 = -(r1+A(1,1))-1/(2*sigma1^2)*(r1-sigma1^2/2)^2;
    alpha2 = -(r2+A(2,2))-1/(2*sigma2^2)*(r2-sigma2^2/2)^2;
    xi = alpha2-alpha1;
    eta = beta2-beta1;
    
    % alpha是时间空间剖分网比，alpha=dt/(dx^2)
    alpha = dt/(dx^2);
    
    % 形成矩阵G
    x = linspace(-L,L,Nx+1);
    g1 = exp(-(alpha1*(t0+(n+1)*dt)+beta1*x)).*max(K-exp(x),0);
    G1 = g1(2:1:Nx).';
    g2 = exp(-(alpha2*(t0+(n+1)*dt)+beta2*x)).*max(K-exp(x),0);
    G2 = g2(2:1:Nx).';
    G = [G1;G2];
    
    % w1,w2是制度1的时间第n+1层的两个端点值，w3,w4是制度2的时间第n+1层的两个端点值
    w1 = exp(-(alpha1*(t0+(n+1)*dt)+beta1*(-L)))*max(K-exp(-L),0);
    w3 = exp(-(alpha2*(t0+(n+1)*dt)+beta2*(-L)))*max(K-exp(-L),0);
    w2 = 0;
    w4 = 0;
    
    % 形成矩阵B1,B2
    B1 = zeros(Nx-1);
    B2 = zeros(Nx-1);
    for i = 2:Nx-2
        B1(i,i-1) = -alpha*sigma1^2/2;
        B1(i,i) = 1+alpha*sigma1^2;
        B1(i,i+1) = -alpha*sigma1^2/2;
        B2(i,i-1) = -alpha*sigma2^2/2;
        B2(i,i) = 1+alpha*sigma2^2;
        B2(i,i+1) = -alpha*sigma2^2/2;
    end
    B1(1,1) = 1+alpha*sigma1^2;
    B1(1,2) = -alpha*sigma1^2/2;
    B1(Nx-1,Nx-2) = -alpha*sigma1^2/2;
    B1(Nx-1,Nx-1) = 1+alpha*sigma1^2;
    B2(1,1) = 1+alpha*sigma2^2;
    B2(1,2) = -alpha*sigma2^2/2;
    B2(Nx-1,Nx-2) = -alpha*sigma2^2/2;
    B2(Nx-1,Nx-1) = 1+alpha*sigma2^2;
    
    % 形成矩阵C1,C2
    a1 = ones(1,Nx-1);
    a2 = ones(1,Nx-1);
    for i = 1:Nx-1
        a1(i) = exp(i*dx*eta);
    end
    for i = 1:Nx-1
        a2(i) = exp(-(i*dx*eta));
    end
    C1 = diag(a1)*A(1,2)*exp(xi*(t0+dt*(n))+eta*x0);
    C2 = diag(a2)*A(2,1)*exp(-(xi*(t0+dt*(n))+eta*x0));
    C1 = dt*C1;
    C2 = dt*C2;
    
    % 形成M和F(全隐格式)
    M = [B1 C1;C2 B2];   
    F = zeros(2*Nx-2,1);
    F(1) = -alpha*sigma1^2/2*w1;     
    F(Nx-1) = -alpha*sigma1^2/2*w2;
    F(Nx) = -alpha*sigma2^2/2*w3;
    F(2*Nx-2) = -alpha*sigma2^2/2*w4;
    F = F-W;
    F = F+M*G;
    
    % 形成M和F(半显半隐，本制度显式，他制度隐式)
%     M = [B1 zeros(size(B1));zeros(size(B1)) B2];
%     D = [-diag(ones(1,Nx-1)) C1;C2 -diag(ones(1,Nx-1))];
%     F = zeros(2*Nx-2,1);
%     F(1) = -alpha*sigma1^2/2*w1;  
%     F(Nx-1) = -alpha*sigma1^2/2*w2;
%     F(Nx) = -alpha*sigma2^2/2*w3;
%     F(2*Nx-2) = -alpha*sigma2^2/2*w4;
%     F = F+D*W;
%     F = F+M*G;
    
end