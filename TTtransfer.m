 function [value1,value2] = TTtransfer(K,Nt,s,T,r1,r2,sigma1,sigma2,A)
    % �ƶ�ת�Ƶ�������ģ�ͣ����������ƶ��·ֱ���ͬһ���㴦��ֵ

    % n��ʱ�����
    % Nt��ʱ��ֻ�����
    % s�ǿռ�����
    % K r1,r2 sigma1,sigma2�Ǹ�����ͬ�ƶ��µĲ���
    % T��ʱ��߽�
    % A���ƶ�ת������

    dt = T/Nt;
    sig = max(sigma1,sigma2)+(sqrt(1.5)-1)*(sigma1+sigma2)/2;
    u = exp(sig*sqrt(dt));
    d = 1/u;
    lambda1 = sig/sigma1;
    lambda2 = sig/sigma2;
    
    % pij�Ľ��ͣ�
    % j=1,2 ��ʾ��һ�������ƶ�
    % i=1,2,3 ��ʾ�۸񱣳����������䡢�½�
    p11 = (exp(r1*dt)-d-(1-1/lambda1^2)*(1-d))/(u-d);
    p12 = (exp(r2*dt)-d-(1-1/lambda2^2)*(1-d))/(u-d);
    p21 = 1-1/lambda1^2;
    p22 = 1-1/lambda2^2;
    p31 = (u-exp(r1*dt)-(1-1/lambda1^2)*(u-1))/(u-d);
    p32 = (u-exp(r2*dt)-(1-1/lambda2^2)*(u-1))/(u-d);

    % �γ�״̬ת�ƾ���P
    I0 = [1,1];
    I = diag(I0);
    P = I+dt*A+dt^2*A^2/2; 
    
    V1 = zeros(Nt+1,2*Nt+1);
    V2 = zeros(Nt+1,2*Nt+1);
    for j = 1:2*Nt+1
        V1(Nt+1,j) = max(K-s*u^Nt*d^(j-1),0);
        V2(Nt+1,j) = max(K-s*u^Nt*d^(j-1),0);
    end
    
    for i = Nt:-1:1
        for j = 1:1:2*i-1
            V1(i,j)=max(K-s*u^(i-1)*d^(j-1),exp(-r1*dt)*(P(1,1)*(p11*V1(i+1,j)+p21*V1(i+1,j+1)+p31*V1(i+1,j+2))+P(1,2)*(p11*V2(i+1,j)+p21*V2(i+1,j+1)+p31*V2(i+1,j+2))));
            V2(i,j)=max(K-s*u^(i-1)*d^(j-1),exp(-r2*dt)*(P(2,1)*(p12*V1(i+1,j)+p22*V1(i+1,j+1)+p32*V1(i+1,j+2))+P(2,2)*(p12*V2(i+1,j)+p22*V2(i+1,j+1)+p32*V2(i+1,j+2))));
        end
    end
    
    value1=V1(1,1);
    value2=V2(1,1);

 end