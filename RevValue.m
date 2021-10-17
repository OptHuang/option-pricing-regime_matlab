function [V1,V2] = RevValue(n,Nx,Phi_end,G,beta1,beta2,alpha1,alpha2,t0,x0,dt,dx)
    % ����Ӻ������ڽ�pcm�㷨��õĵ�n+1��Ľ⻹ԭ������Ҫ�ĺ���V����Ӧ�����ֵ
    
    % n��ָ��ʱ��ֻ���n�㣬Nx�ǿռ�ֻ�����
    % Phi_end��pcmֱ����õ�(2Nx-2)*1��������
    % G�Ǹò��Ӧ��(2Nx-2)*1��������
    % V1,V2�����ƶ�1,2�·ֱ��Ӧ��Value��������ȡ�ĵ㴦��ֵ
    % beta1,beta2,alpha1,alpha2��AssembleMatrix���������ת���еĲ���
    % t0,x0,dt,dx�ֱ���ʱ����㣬�ռ���㣬ʱ��ֻ�ÿ�γ��ȣ��ռ�ֻ�ÿ�γ���
    
    Phi = Phi_end+G;
    W1 = Phi(1:Nx-1);
    W2 = Phi(Nx:2*Nx-2);
    
    c1 = zeros(1,Nx-1);
    c2 = zeros(1,Nx-1);
    for i =1:Nx-1
        c1(i) = exp(alpha1*(t0+(n+1)*dt)+beta1*(x0+i*dx));
        c2(i) = exp(alpha2*(t0+(n+1)*dt)+beta2*(x0+i*dx));
    end
    reC1 = diag(c1);
    reC2 = diag(c2);
    
    V1 = reC1*W1;
    V2 = reC2*W2;

end