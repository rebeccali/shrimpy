function [Fb, Mb, aoa1] = blade(Vcg_f,omgBlade_f,nu,pitch1,Cl, Cd, h, beta, chord, span, dSpan, rho)
    % vcg in the flyer frame
    %omgBlade_f is the omg in flyer frame plus the blade omega in the
    %frame
    % nu is inflow velocity
    % h is height of prop to cg??
    % span is a vector from 0 to RMax incremented by dSpan 
    drSteps = length(beta);

    Fb = zeros(drSteps,3);
    Mb = zeros(drSteps,3);

    %calculate relative wind
    Ut1 = Vcg_f(1) + omgBlade_f(2)*h - span*omgBlade_f(3); %correct
    Up1 = Vcg_f(3) + span*omgBlade_f(1) - nu;
    rw = atan2(Up1,Ut1);

    %lift and drag
    aoa1=mod(pi+beta+pitch1+rw,2*pi)-pi;
    ind = int32((aoa1+pi)*1000000+1);
    qs = rho/2*(Ut1.*Ut1 + Up1.*Up1).*chord;
    l1 = qs.*Cl(ind)';
    d1 = qs.*abs(Cd(ind))'; 

    %p and t components of force (n is up, t is back)
    n1 = l1.*cos(rw)+d1.*sin(rw);
    t1 = -l1.*sin(rw)+d1.*cos(rw);

    %TODO::should be done for between elements, not each element
    %force vector
    Fb(:,1) = -t1*dSpan;
    Fb(:,3) = -n1*dSpan;

    %moment vector
    Mb(:,1) = -span.*n1*dSpan;
    Mb(:,2) = -h.*t1*dSpan; %correct (sign change for h)
    Mb(:,3) = span.*t1*dSpan;
end