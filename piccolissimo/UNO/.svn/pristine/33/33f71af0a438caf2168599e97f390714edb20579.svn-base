function nu_n = NumericalFixedPointInflow(nu_n)

global Vcg T rho span_d dSpan_d nu

    alpha_base = .1;
    V_inf = sqrt(Vcg(1)*Vcg(1) + Vcg(2)*Vcg(2));
    nu = inf;
    iterations = 0;
    
%     while(nnz(abs((nu_n-nu)./nu_n) > .0005*alpha_base))
    while(nnz(abs((nu_n-nu)) > .1*alpha_base) && iterations < 1000)
        nu = nu_n;
        [F_p_nu, ~, F_d_nu, ~] = ComputeAero();
        T = -F_p_nu(:,3)-F_d_nu(:,3);
        nu_t = sqrt(T'./(2*rho*pi*(span_d.*span_d-(span_d-dSpan_d).*(span_d-dSpan_d)))); % Assumes span_d = span_p!
        alpha = alpha_base*(1-exp(-abs(nu_t))).^5;
        nu_n = nu_n.*(1-alpha) + alpha.*(nu_t.*nu_t)./sqrt(V_inf*V_inf+nu_n.*nu_n);
        nu_n(isnan(nu_n)|T'==0) = 0;
        iterations = iterations + 1;
    end
    nu = nu_n;
    
end