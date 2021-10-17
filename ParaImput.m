function [T,L,Nt,Nx,dx,dt,sigma1,sigma2,r1,r2,gamma,mu,eps,rho0,A,x0,t0,K] = ParaImput()
    % 由于参数比较多，所以写一个传入参数的函数方便修改调用
    
    T = 1;  % 传入时间上界T
    K = 9;  % 传入K值 
    
    % 传入两种制度下的波动率sigma和无风险利率r
    sigma1 = 0.8;
    sigma2 = 0.3;
    r1 = 0.1;
    r2 = 0.05;
    
    % 传入pcm的几个参数
    gamma = 0.9;
    mu = 0.1;
    rho0 = 1.9;
    eps = 1e-6;
    
    % 截取边界
    L = XspaceBoundry(sigma1,sigma2,r1,r2,T,K); % 计算出空间边界L
    
    % 传入不同制度之间的转化矩阵A
    a11 = 6;
    a12 = -6;
    a21 = -9;
    a22 = 9;
    A = [a11 a12;a21 a22];
    
    x0 = -L;  % 传入空间左侧起点
    t0 = 0;   % 传入时间起点
    
    Nt = 300;  % 传入时间分划段数
    Nx = 150;   % 传入空间分划段数

    dt = T/Nt;    % 计算出时间分划间距
    dx = 2*L/Nx;  % 计算出空间分划间距
    
end