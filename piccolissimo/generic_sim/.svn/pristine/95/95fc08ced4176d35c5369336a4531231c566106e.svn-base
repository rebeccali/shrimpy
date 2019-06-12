function [F1 M1] = blade_mod(Vcg,omg,nu,pitch1,bladeProperties,rho)
    h = bladeProperties(1);
    R = bladeProperties(2);
    drSteps = bladeProperties(3);
    beta = bladeProperties(4:3+drSteps);
    chord = bladeProperties(4+drSteps:3+drSteps+drSteps);
    % persistent Cd Cl
    % if isempty(Cd)
    %     load('Cd');
    %     load('Cl');
    % end
    persistent Cd_segmented Cl_segmented
    if isempty(Cd_segmented)
        load('Cd_segmented');
        load('Cl_segmented');
    end

    FM = [0 0 0 0 0 0];
    options = odeset('RelTol',1e-1, 'AbsTol', 1e-3);
    [~, FMout] = ode113(@bladeODE,linspace(0,R,drSteps),FM,options);
    
    function dFM = bladeODE(r,~)
        %calculate relative wind
        Ut1 = Vcg(1) - (omg(2)*h+r*omg(3));
        Up1 = Vcg(3) + r*omg(1) - nu;

        %lift and drag
%         betaCurrent = interp1q(linspace(0,R,drSteps),beta',r);
        betaCurrent = beta(int32(r/R*(drSteps-1) + 1));
%         chordCurrent = interp1q(linspace(0,R,drSteps),chord',r);
        chordCurrent = chord(int32(r/R*(drSteps-1) + 1));
        aoa1=mod(pi+betaCurrent+pitch1+atan2(Up1,Ut1),2*pi)-pi;
        l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chordCurrent*abs(Cl_segmented(int32((aoa1+pi)*100+1)));
        d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chordCurrent*abs(Cd_segmented(int32((aoa1+pi)*100+1)));
    %     l1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(interp1q(Cl(:,1),Cl(:,2),aoa1));
    %     d1 = rho/2*(Ut1*Ut1 + Up1*Up1)*chord(i)*abs(interp1q(Cd(:,1),Cd(:,2),aoa1));

        %p and t components of force (n is up, t is back)
        n1 = l1*cos(atan2(Up1,Ut1))+d1*sin(atan2(Up1,Ut1));
        t1 = -l1*sin(atan2(Up1,Ut1))+d1*cos(atan2(Up1,Ut1));

        dFM = [-t1; 0; -n1; -r*n1; h*t1; r*t1];
    end
    F1(1) = FMout(end,1);
    F1(2) = FMout(end,2);
    F1(3) = FMout(end,3);
    M1(1) = FMout(end,4);
    M1(2) = FMout(end,5);
    M1(3) = FMout(end,6);
end