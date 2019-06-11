figure(1)
scaleFactorZ = 1e-3;
scaleFactorY = 17e-3;
scaleFactorX = 4e-3;
[F,V] = stlread('file.STL');
bodySFZ = 6e-3;
bodySFY = 2e-3;
bodySFX = 2e-3;

bound = 50e-3;
mass = 1;
%make ground plane
groundX = [-10e3, 10e3, 10e3, -10e3];
groundY = [-10e3, -10e3, 10e3, 10e3];
%propeller vertices definiton
x = [0,1,1,0,0,0;1,1,0,0,1,1;1,1,0,0,1,1;0,1,1,0,0,0];
y = [0,0,1,1,0,0;0,1,1,0,0,0;0,1,1,0,1,1;0,0,1,1,1,1];
z = [0,0,0,0,0,1;0,0,0,0,0,1;1,1,1,1,0,1;1,1,1,1,0,1];

%body vertices definiton
bodyX = [0,1,1,0,0,0;1,1,0,0,1,1;1,1,0,0,1,1;0,1,1,0,0,0];
bodyY = [0,0,1,1,0,0;0,1,1,0,0,0;0,1,1,0,1,1;0,0,1,1,1,1];
bodyZ = [0,0,0,0,0,1;0,0,0,0,0,1;1,1,1,1,0,1;1,1,1,1,0,1];

a = 0.5*ones(4,6);
x = x-a;
y = y-a;
z = z-a; 

bodyX = bodyX-a;
bodyY = bodyY-a;
bodyZ = bodyZ-2*a;

z = z.*scaleFactorZ;
y = y.*scaleFactorY;
x = x.*scaleFactorX;

bodyZ = bodyZ.*bodySFZ;
bodyY = bodyY.*bodySFY;
bodyX = bodyX.*bodySFX;

ground = patch(groundX, groundY, 'red');
body = patch(bodyX,bodyY,bodyZ, 'green');
alpha(0.3);
vertexGround = ground.Vertices;
myProp = patch(x,y,z,'blue');
%myProp = patch(V(:,1)./1000,V(:,2)./1000,V(:,3)./1000,'blue');
vertex = myProp.Vertices;
vertexOrig = myProp.Vertices;
vertexBody = body.Vertices;
vertexBodyOrig = body.Vertices; 
grid on;

xlim([-bound,bound]);
ylim([-bound,bound]);
zlim([-20e-3, 3*bound]);





%timestep (s)
dt = 0.1/6;
d = [0,0,0]';
transX = 0;
transY = 0;
transZ = 0;

vertexNew = ones(24,4);
vertexBodyNew = ones(24,4);

[t,states] = ode45(@myfun,[0 20], [0;0;0;0.8;0;0]);


tAnim = (0:dt:20)';
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

%% moment calculation stuff

%force and distance are 3x1 vectors 
% t = fxl 



%% animation test
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










