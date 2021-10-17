function [M,F,G,beta1,beta2,alpha1,alpha2] = AssembleMatrix(n,Nx,dx,dt,sigma1,sigma2,r1,r2,A,x0,t0,W,L,K)
    % �γɾ���M��F

    % n���ֽ׶ε����������ڵ�ʱ�䷽��Ĳ���(n�㵽n+1��)
    % Nx�ǿռ䷽���ʷֶ�����dx��ÿС�γ��ȣ�x0��������
    % dt��ʱ��ֻ�ÿС�γ��ȣ�t0�����·���ʼ��
    % sigma1,sigma2,r1,r2�Ǹ���ϵ��,A�Ǹ�����ת������
    % W��һ��(2Nx-2)*1������������ʱ����n�������
    % L�ǽ�ȡ�߽磬K�ǳ�����
    
    % ����任���µĲ����仯
    beta1 = 1/2-r1/(sigma1^2);
    beta2 = 1/2-r2/(sigma2^2);
    alpha1 = -(r1+A(1,1))-1/(2*sigma1^2)*(r1-sigma1^2/2)^2;
    alpha2 = -(r2+A(2,2))-1/(2*sigma2^2)*(r2-sigma2^2/2)^2;
    xi = alpha2-alpha1;
    eta = beta2-beta1;
    
    % alpha��ʱ��ռ��ʷ����ȣ�alpha=dt/(dx^2)
    alpha = dt/(dx^2);
    
    % �γɾ���G
    x = linspace(-L,L,Nx+1);
    g1 = exp(-(alpha1*(t0+(n+1)*dt)+beta1*x)).*max(K-exp(x),0);
    G1 = g1(2:1:Nx).';
    g2 = exp(-(alpha2*(t0+(n+1)*dt)+beta2*x)).*max(K-exp(x),0);
    G2 = g2(2:1:Nx).';
    G = [G1;G2];
    
    % w1,w2���ƶ�1��ʱ���n+1��������˵�ֵ��w3,w4���ƶ�2��ʱ���n+1��������˵�ֵ
    w1 = exp(-(alpha1*(t0+(n+1)*dt)+beta1*(-L)))*max(K-exp(-L),0);
    w3 = exp(-(alpha2*(t0+(n+1)*dt)+beta2*(-L)))*max(K-exp(-L),0);
    w2 = 0;
    w4 = 0;
    
    % �γɾ���B1,B2
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
    
    % �γɾ���C1,C2
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
    
    % �γ�M��F(ȫ����ʽ)
    M = [B1 C1;C2 B2];   
    F = zeros(2*Nx-2,1);
    F(1) = -alpha*sigma1^2/2*w1;     
    F(Nx-1) = -alpha*sigma1^2/2*w2;
    F(Nx) = -alpha*sigma2^2/2*w3;
    F(2*Nx-2) = -alpha*sigma2^2/2*w4;
    F = F-W;
    F = F+M*G;
    
    % �γ�M��F(���԰��������ƶ���ʽ�����ƶ���ʽ)
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