function L = XspaceBoundry(sigma1,sigma2,r1,r2,T,K)
    % �������ضϱ߽糤�ȵĺ���

    % L0,X0����һ��2*1��������
    L0 = zeros(2,1);
    L0(1) = -1.25*sigma1^2*T*(r1/sigma1^2-0.5)+0.5*sqrt(25/4*sigma1^4*T^2*(r1/sigma1^2-0.5)^2-10*sigma1^2*T*log(1e-6/sqrt(5*K)));
    L0(2) = -1.25*sigma2^2*T*(r2/sigma2^2-0.5)+0.5*sqrt(25/4*sigma2^4*T^2*(r2/sigma2^2-0.5)^2-10*sigma2^2*T*log(1e-6/sqrt(5*K)));
    
    X0 = zeros(2,1);
    
    r = [r1 r2];
    sigma = [sigma1 sigma2];
    
    % ��߽�X_0
    for i =1:2
        X0(i) = 2*r(i)/(2*r(i)+sigma(i)^2);
    end
    if (X0(1)<X0(2))
        X_0 = X0(1);
    else
        X_0 = X0(2);
    end
    
    % �ұ߽�L_0
    if (L0(1)>L0(2))
        L_0 = L0(1);
    else
        L_0 = L0(2);
    end
    
    % ���ձ߽�L
    if (-log(X_0)>L_0)
        L = -log(X_0);
    else
        L = L_0;
    end
    
end