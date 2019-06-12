function [tout Xout] = asianTop(X,time)
    %Syntax:
    %X = [nu; u; v; w; p; q; omg; phi; theta; psi; xe; ye; ze];
    %Example:
    %X = [3; 0; 0; 0; 0; 0; -200; .1; 0; 0; 0; 0; 0];
    %asianTop(X,3);

    %Do math
    [tout Xout] = ode113(@asianTopODE,[0 time],X);
    
    %Make plots
    figure(1);
    hold off;
    plot(tout,Xout(:,11),'color','r'); %x
    hold on;
    title('Position');
    plot(tout,Xout(:,12),'color','g'); %x
    plot(tout,Xout(:,13),'color','b'); %x
    figure(2); 
    hold off;
    plot(tout,Xout(:,8),'color','r'); %angle rot
    hold on;
    title('Angle');
    plot(tout,Xout(:,9),'color','g'); %angle rot
%     plot(tout,Xout(:,10),'color','b'); %angle rot
    figure(3);
    hold off;
    plot(tout,Xout(:,2),'color','r');%velocity up
    hold on;
    title('Velocity');
    plot(tout,Xout(:,3),'color','g');%velocity up
    plot(tout,Xout(:,4),'color','b');%velocity up
    figure(4);
    hold off;
    plot(tout,Xout(:,5),'color','r');%rotation
    hold on;
    title('Angular Velocity');
    plot(tout,Xout(:,6),'color','g');%rotation
    plot(tout,Xout(:,7),'color','b');%rotation
    figure(5);
    hold off;
    plot(tout,Xout(:,1)); %nu
    hold on;
    title('Nu');
    sz = size(tout);
    orientation = zeros(sz(1),3);
    for i = 1:sz(1)
        angles(1) = Xout(i,8);
        angles(2) = Xout(i,9);
        angles(3) = Xout(i,10);
        Leb = [cos(angles(2))*cos(angles(3)),sin(angles(1))*sin(angles(2))*cos(angles(3))-cos(angles(1))*sin(angles(3)),cos(angles(1))*sin(angles(2))*cos(angles(3))+sin(angles(1))*sin(angles(3)); cos(angles(2))*sin(angles(3)),sin(angles(1))*sin(angles(2))*sin(angles(3))+cos(angles(1))*cos(angles(3)),cos(angles(1))*sin(angles(2))*sin(angles(3))-sin(angles(1))*cos(angles(3));-sin(angles(2)),sin(angles(1))*cos(angles(2)),cos(angles(1))*cos(angles(2))];
        orientation(i,:) = (Leb*[0 0 1]')';
    end
    figure(6);
    hold off;
    plot(tout,orientation(:,1),'color','r');%rotation
    hold on;
    title('orientation');
    plot(tout,orientation(:,2),'color','g');%rotation
    plot(tout,orientation(:,3),'color','b');%rotation
end


function dX = asianTopODE(t,X)
    
persistent h_p R drSteps beta chord Dp mProp pArea H bladeProperties Db mBody m rho g tauI kp ki kd deltaPrev iClamp tau I

if isempty(h_p)
    %Prop properties
    h_p = .03; %Height of blade over cg (cg assumed axial) (meters)
    R = .1; %Blade radius (meters)
    drSteps = 10; %number of changes in beta and chord wrt r (1 = constant)
    beta = 0.5*ones(1,drSteps); %Blade twist (relative to zero lift)(assumed constant) (radians)
    %beta = [? ? ? ?] to be function of r
    chord = .02*ones(1,drSteps); %Chord lengh (assumed constant) (meters)
    %chord =[? ? ? ?] to be function of r

    Dp = [.02,.2,.001]; %xyz dimentions of the blade for interta
    mProp = .005; %prop mass in kg
    pArea = R*R*pi;
    H = R; %approx height above prop that air is moved (estimated to be radius of prop)
    bladeProperties(1) = h_p;
    bladeProperties(2) = R;
    bladeProperties(3) = drSteps;
    bladeProperties = [bladeProperties beta];
    bladeProperties = [bladeProperties chord];
    
    %Body properties
    h_b = -.056;
    Db = [.005,.15]; %r,L dimentions of the body
    mBody = .006; %body mass in kg
    
    m = mProp+mBody;
    
    %Environment
    rho = 1.225; %air density (sea level 1.225) (kg/m^3)
    g = 9.8; %gravity acceleration (m/s^2)
    
    %Power setup
%     tauI = .0013;
%     kp = .001;
%     ki = .001;
%     kd = 5;
%     deltaPrev = 0;
%     % iClamp = .005; %max power
%     iClamp = 0; %unpowered flight
    tau = 0;

    Ibs = [1/12*mBody*(3*Db(1)*Db(1)+Db(2)*Db(2)),0,0;0,1/12*mBody*(3*Db(1)*Db(1)+Db(2)*Db(2)),0;0,0,mBody*Db(1)*Db(1)/2];
    Iob = mBody*(h_b*h_b*eye(3)+[0; 0; h_b]*[0, 0, h_b]);
    Ips = [1/12*mProp*(Dp(2)*Dp(2)+Dp(3)*Dp(3)),0,0;0,1/12*mProp*(Dp(1)*Dp(1)+Dp(3)*Dp(3)),0;0,0,1/12*mProp*(Dp(1)*Dp(1)+Dp(2)*Dp(2))];
    Iop = mProp*(h_p*h_p*eye(3)+[0; 0; h_p]*[0, 0, h_p]);
    I = Ips+Iop+Ibs+Iob;

end

    nu = X(1);
    Vcg = X(2:4);
    omg = X(5:7);
    angles = X(8:10);
    xcg = X(11:13);

    pitch1 = 0;
    pitch2 = 0;

    %Thrust
%     [F1 M1 F2 M2] = blades(Vcg,omg,nu,pitch1,pitch2,bladeProperties,rho);
    Rr_f1 = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    [F1 M1] = blade(Rr_f1*Vcg,Rr_f1*omg,nu,pitch1,bladeProperties,rho);
    F1 = (F1*Rr_f1);
    M1 = (M1*Rr_f1);
    
    Rr_f2 = [-1, 0, 0; 0, -1, 0; 0, 0, 1];
    [F2 M2] = blade(Rr_f2*Vcg,Rr_f2*omg,nu,pitch2,bladeProperties,rho);
    F2 = (F2*Rr_f2);
    M2 = (M2*Rr_f2);
    
    T = -F1(3)-F2(3);
    
%     %compute desired torque
%     delta = xcg(3)+5;
%     tauI = tauI + ki*delta*dt;
%     tau = tauI+delta*kp + (delta-deltaPrev)*kd;
%     deltaPrev = delta;
%     if tauI < 0
%         tauI = 0;
%     end
%     if tauI > iClamp
%         tauI = iClamp;
%     end
%     if tau < 0
%         tau = 0;
%     end
%     if tau > iClamp
%         tau = iClamp;
%     end

    %Change in induced air velocity
    nd1 = -2*nu*abs(-Vcg(3) + nu)/H;
    nd2 = T/(rho*pArea*H);
    nuDot = nd1+nd2;
    %Force caused by gravity
    Fg = [-m*g*sin(angles(2)), m*g*cos(angles(2))*sin(angles(1)), m*g*cos(angles(2))*cos(angles(1))];
    %CG acceleration
    VcgDot = (1/(mProp+mBody))*(F1+F2+Fg)'-cross(omg,Vcg);
    %CG rotational acceleration
    omgDot = I\(M1+M2+[0,0,-tau]-cross(omg',I*omg))';
    %World rotational acceleration
    anglesDot = [omg(1)+(omg(2)*sin(angles(1))+omg(3)*cos(angles(1)))*tan(angles(2)); omg(2)*cos(angles(1))-omg(3)*sin(angles(1));(omg(2)*sin(angles(1))+omg(3)*cos(angles(1)))*sec(angles(2))];
    %World velocity
    Leb = [cos(angles(2))*cos(angles(3)),sin(angles(1))*sin(angles(2))*cos(angles(3))-cos(angles(1))*sin(angles(3)),cos(angles(1))*sin(angles(2))*cos(angles(3))+sin(angles(1))*sin(angles(3)); cos(angles(2))*sin(angles(3)),sin(angles(1))*sin(angles(2))*sin(angles(3))+cos(angles(1))*cos(angles(3)),cos(angles(1))*sin(angles(2))*sin(angles(3))-sin(angles(1))*cos(angles(3));-sin(angles(2)),sin(angles(1))*cos(angles(2)),cos(angles(1))*cos(angles(2))];
    XeDot = Leb*Vcg;
    
    dX = [nuDot; VcgDot(1); VcgDot(2); VcgDot(3); omgDot(1); omgDot(2); omgDot(3); anglesDot(1); anglesDot(2); anglesDot(3); XeDot(1); XeDot(2); XeDot(3)];
end