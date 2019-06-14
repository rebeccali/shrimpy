function [ Xret, found_trim] = FlyerTrim(time)
    global Xbase m_s m_r R_nu rho
    
    warnId = 'MATLAB:ode45:IntegrationTolNotMet';
    warning('error', warnId);    
    
    Xbase(1) = MomentumInflow((m_s+m_r), R_nu, rho);
    args = {'throttle', 'trim','plot'};
%     try     
    [tout, Xout] = FlyerSimulation(Xbase,time,args);

    Xret = Xout(end,:)';
    Xret(13) = 0;
    Xret(14) = 0;
    if(tout(end) == time)
        found_trim = false;
    else
        found_trim = true;
    end
    
%     catch trim_error
%         % Check if we indeed failed on meeting the tolerances
%         if strcmp(trim_error.identifier, warnId)
%             found_trim = false;
%         else
%             % Something else has gone wrong: just re-throw the error
%             throw(trim_error);
%         end
% 
%     end
end
