function [ Xret ] = findTrim(time)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    global waypoints g Xbase
    
%     if isempty(g)
%         setup_flyer_ideal_prop_on_bottom;
%     end
    
    waypoints = [0,0,0];
    if (exist('Xbase','var') && ~isempty(Xbase))
        X = Xbase
    else
    %     X = [3.1426; 0; 0; 0; 0; 0; 0; -297.5032; 37.6379; 0; 0; 0; 0; 0; 0; 0; 0; -2.3980];
%         X = [4.8883; 0; 0; 0; 0; 0; 0; -549.4082; 34.6068; 0; 0; 0; 0; 0; 0; 0; 0; -2.4730]
    %     X = [4.8883; 0; 0; 0; 0; 0; 0; -3.628421450098252e+02; 31.106515996471420; 0; 0; 0; 0; 0; 0; 0; 0; -1.027565017];
        X = [2.765; 0; 0; 0; 0; 0; 0; 502.8971; -28.5059; 0; 0; 0; 0; 0; 0; 0; 0; 5.6877]
    end
    args = {'throttle', 'trim'};
    [tout Xout] = yimFlyerLite(X,time,args);
    Xret = Xout(end,:)';
    Xret(13) = 0;
    Xret(14) = 0;
end

