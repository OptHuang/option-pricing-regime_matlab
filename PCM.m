function Phi_end = PCM(Nx,eps,gamma,mu,rho0,M,F,Phi_start)
    % pcm�㷨���ĳ���

    % Nx�ǿռ䷽���ʷֶ���
    % eps����ֹ�������
    % gamma��mu��rho0�Ǹ�����(0,1)֮���ϵ��
    % M�Ǹ�����(2Nx-2)*(2Nx-2)��ϵ������
    % F�Ǹ�����(2Nx-2)*1��������
    % Phi_start�Ǹò������ĳ�ʼֵ
    
    Phi = Phi_start;
    
    beta = 1;
    k = 0;
    F_u = M*Phi+F;
    tol = norm(Phi-max(Phi-F_u,zeros(2*Nx-2,1)),inf);
    while (tol>eps)
        % ͶӰ
        Phi_mid = Phi;
        F_u_mid = F_u;
        Phi = max(Phi_mid-beta.*F_u_mid,0);
        % ���з���
        F_u = M*Phi+F;
        du = Phi_mid-Phi;
        dF = beta.*(F_u_mid-F_u);
        rho = norm(dF)/norm(du);
        while (rho>gamma)
            beta = 2/3*beta*min(1,1/rho);
            Phi = max(Phi_mid-beta.*F_u_mid,zeros(2*Nx-2,1));
            F_u = M*Phi+F;
            du = Phi_mid-Phi;
            dF = beta*(F_u_mid-F_u);
            rho = norm(dF)/norm(du);
        end
        duF = du-dF;
        r1 = dot(du,duF);
        r2 = dot(duF,duF);
        alpha = r1/r2;
        Phi = Phi_mid-alpha*rho0*duF;
        F_u = M*Phi+F;
        tol = norm(Phi-max(Phi-F_u,zeros(2*Nx-2,1)),inf);
        if (rho<mu)
           beta = beta*rho0; 
        end
        k = k+1;
    end     
    Phi_end = Phi;
    
end