function Pic3DPCM
    % 绘制三维曲面(PCM)
    
    format long
    
    [T,L,Nt,Nx,dx,dt,sigma1,sigma2,r1,r2,gamma,mu,eps,rho0,A,x0,t0,K] = ParaImput();
    [Utrue1,Utrue2] = DataMatrix(Nt,Nx,dx,dt,sigma1,sigma2,r1,r2,gamma,mu,eps,rho0,A,x0,t0,K);
    [S,T] = MeshGeneration(T,L,Nx,Nt);
    
    % 现实中矩阵元素对应坐标需要将矩阵以空间方向为轴翻转一下
    Utrue1 = flipud(Utrue1);
    Utrue2 = flipud(Utrue2);
    
    subplot(1,3,1)
    mesh(S,T,Utrue1)
    axis([0 25 0 1 0 10])
    
    subplot(1,3,2)
    mesh(S,T,Utrue2)
    axis([0 25 0 1 0 10])
    
    subplot(1,3,3)
    x = linspace(-L,L,Nx+1);
    s = exp(x);
    plot(s,Utrue1(205,:))
    axis([0 25 0 9])
    
end