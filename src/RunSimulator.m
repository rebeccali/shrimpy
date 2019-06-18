%% Runs a simulation of the vehicle and body dynamics 

%% Clean Space
close all;
clear all;

%% Define Simulation Parameters
dt = 0.01;
startTime = 0;
endTime = 20;
x0 = [0;0;0;0.8;0;0];

%% Do Simulation

[t,states] = ode45(@FlyerDynamics, [startTime, endTime], x0);

%% Set Up figures 
figure;
grid on;
bounds = 50e-3;
xlim([-bounds,bounds]);
ylim([-bounds,bounds]);
zlim([-20e-3, 3*bounds]);


% make ground plane
groundX = [-10e3, 10e3, 10e3, -10e3];
groundY = [-10e3, -10e3, 10e3, 10e3];
ground = patch(groundX, groundY, 'red');

% Set up body shape
flyerVerticies = getFlyerBodyVerticies();
body = patch(flyerVerticies.bodyX, flyerVerticies.bodyY, flyerVerticies.bodyZ, 'green');
myProp = patch(flyerVerticies.propX, flyerVerticies.propY, flyerVerticies.propZ,'blue');

alpha(0.3); % make ground transparent
vertexGround = ground.Vertices;
vertex = myProp.Vertices;
vertexOrig = myProp.Vertices;
vertexBody = body.Vertices;
vertexBodyOrig = body.Vertices; 


%% Do Animation
tAnim = (startTime:dt:endTime)';
thetaAnim = interp1(t, states(:,1),tAnim);
thetaAnim2 = interp1(t, states(:,3),tAnim);
thetaAnim3 = interp1(t, states(:,5),tAnim);

figure(2);
plot(t,states(:,1),'b');
hold on;
plot(t,states(:,2),'r');
legend('theta (rad)', 'omega (rad/s)');

%establish unit vector for prop thrust direction 
%initialized as [0 0 1]
propX = 0;
propY = 0; 
propZ = 1;
propDirec = [propX propY propZ];
%motor torque input plot
%represent as a square wave

transX = 0;
transY = 0;
transZ = 0;

vertexNew = ones(24, 4);
vertexBodyNew = ones(24, 4);

d = [0,0,0]';

for i=1:length(tAnim)

    % calculate angles for the orientation of the thrust vector before
    % applying quaternion rotation in the updated axis of rotation 
    % prop always rotates about center 
    % thrust vector along the axis of rotation, normal to the propeller
    % surface 
    
    %use unit vector as the direction for the thrust calculations 
    
    %use new alpha in quaternion to rotate prop
    %use rotation matrix to set the angle for axis of rotation used for quaternion below
    %angle for rotation is determined by the pitching induced from
    %eccentric thrust 
    axisRotateVec = [(cos(thetaAnim3(i)/2)), sin(thetaAnim3(i)/2)*[0,1,0]];
    axisRotateVec = axisRotateVec./sqrt((axisRotateVec(1))^2+(axisRotateVec(2))^2+(axisRotateVec(3))^2);
    a1 = quaternion(axisRotateVec);

    
    qVec = [(cos(thetaAnim(i)/2)), sin(thetaAnim(i)/2)*[0 0 1]];
    
    %normalize quaternion 
    qVec = qVec./(sqrt(qVec(1)^2+qVec(2)^2+qVec(3)^2+qVec(4)^2));
    q1 = quaternion(qVec);
    %convert quaternion to rotation matrix
    rotm = quat2rotm(q1);
    rotm2 = quat2rotm(a1);
    
    H = [rotm, d; 0, 0, 0, 1];
    H2 = [rotm2, [0;0;0]; 0,0,0,1];

    for j = 1:24 
        vertexBodyNew(j,:) = H*H2*[vertexBodyOrig(j,:),1]';
        vertexNew(j,:) = H*H2*[vertexOrig(j,:),1]';
    end

    vertex = vertexNew(:,1:3);
    vertexBody = vertexBodyNew(:,1:3);
    d = [0, 0, thetaAnim2(i)]';
    %apply quaternion
    %vertex = rotatepoint(q1, vertexOrig);
    %draw figure
    figure(1)
    view(3)
    myProp.Vertices =  vertex;
    body.Vertices = vertexBody;
   
    %pause length controls timestep
    pause(dt)
    t = t+dt;
    
    
end

%%










