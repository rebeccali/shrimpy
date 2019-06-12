function [ Xret, found_trim] = FlyerTrim(time)
    global Xbase m_s m_r R_nu rho

    Xbase(1) = MomentumInflow((m_s+m_r), R_nu, rho);
    
    args = {'throttle', 'trim','simple_nu'};
    [tout, Xout] = FlyerSimulation(Xbase,time,args);
    Xret = Xout(end,:)';
    Xret(13) = 0;
    Xret(14) = 0;
    if(tout(end) == time)
        found_trim = false;
    else
        found_trim = true;
    end
end
