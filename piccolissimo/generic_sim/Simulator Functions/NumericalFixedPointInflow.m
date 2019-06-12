function nu_n = NumericalFixedPointInflow(nu_n)

global Vcg nu T rho R_nu

    alpha = .1;
    V_inf = sqrt(Vcg(1)*Vcg(1) + Vcg(2)*Vcg(2));
    flag = true;
%             nu_h = sqrt(m*g/(rho*area_p*2));

    while(abs((nu_n-nu)/nu_n) >= .001 || flag)
        flag = false;
        nu = nu_n;
        [F_p_nu, ~, F_d_nu, ~] = ComputeAero();
        T = -F_p_nu(3)-F_d_nu(3);
        nu_t = sqrt(T/(2*rho*pi*R_nu*R_nu));
        nu_n = V_inf + nu_n*(1-alpha) + alpha*(nu_t*nu_t)/sqrt(V_inf*V_inf+nu_n*nu_n);
        alpha = alpha*.999;
    end

%     if(abs((nu_n1-nu_n)/nu_n1) >= .0005)
%         nu_n1 = NumericalFixedPointInflow(nu_n1);
%     end
end