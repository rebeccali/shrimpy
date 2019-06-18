function [t, states] = animateFlyer(t, states, dt, startTime, endTime)
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
bodyPatch = patch(flyerVerticies.bodyX, flyerVerticies.bodyY, flyerVerticies.bodyZ, 'green');
alpha(0.3); % make ground transparent

propPatch = patch(flyerVerticies.propX, flyerVerticies.propY, flyerVerticies.propZ,'blue');

%% Do Animation
tAnim = (startTime:dt:endTime)';

zetaAnim = interp1(t, states(:,1),tAnim);
zetaAnim2 = interp1(t, states(:,3),tAnim);
zetaAnim3 = interp1(t, states(:,5),tAnim);
zetaAnim4 = interp1(t, states(:,7), tAnim);
zetaAnim5 = interp1(t, states(:,9), tAnim);
zetaAnim6 = interp1(t, states(:,11), tAnim);


%establish unit vector for prop thrust direction 
%initialized as [0 0 1]
propX = 0;
propY = 0; 
propZ = 1;
propDirec = [propX propY propZ];
%motor torque input plot
%represent as a square wave


d = [0,0,0]';
transX = 0;
transY = 0;
transZ = 0;


%% animation
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
    axisRotateVec = [(cos(zetaAnim3(i)/2)), sin(zetaAnim3(i)/2)*[0,0,1]];
    axisRotateVec = axisRotateVec./sqrt((axisRotateVec(1))^2+(axisRotateVec(2))^2+(axisRotateVec(3))^2);
    a1 = quaternion(axisRotateVec);
    
    pitchRotate = [(cos(zetaAnim4(i)/2)), sin(zetaAnim4(i)/2)*[0,1,0]];
    pitchRotate = pitchRotate./sqrt((pitchRotate(1))^2+(pitchRotate(2))^2+(pitchRotate(3))^2);
    p1 = quaternion(pitchRotate);
    
    qVec = [(cos(zetaAnim(i)/2)), sin(zetaAnim(i)/2)*[0 0 1]];
    
    % normalize quaternion 
    qVec = qVec./(sqrt(qVec(1)^2+qVec(2)^2+qVec(3)^2+qVec(4)^2));
    q1 = quaternion(qVec);
    %convert quaternion to rotation matrix
    % TODO: Replace with inbuilt function quat2rotm
    rotm = quat2rotm(q1);
    rotm2 = quat2rotm(a1);
    rotm3 = quat2rotm(p1);
    H = [rotm, d; 0, 0, 0, 1];
    H2 = [rotm2, d; 0,0,0,1];
    H3 = [rotm3, [0;0;0]; 0, 0, 0, 1];
    for j = 1:24 
        vertexBodyNew(j,:) = (H2*H3*[bodyPatch.Vertices(j,:),1]');
        vertexNew(j,:) = H*H3*[propPatch.Vertices(j,:),1]';
    end

    vertex = vertexNew(:,1:3);
    vertexBody = vertexBodyNew(:,1:3);
    d = [zetaAnim5(i), zetaAnim6(i), zetaAnim2(i)]';
    %apply quaternion
    %vertex = rotatepoint(q1, vertexOrig);
    %draw figure
    figure(1)
    view(3) 
    propPatch.Vertices =  vertex;
    bodyPatch.Vertices = vertexBody;
   
    %pause length controls timestep
    pause(dt)
    % TODO: replace this t with something else,
    t = t+dt;
end
end