function animateFlyer(t, states, dt, startTime, endTime)
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

posWorld2Body = [0,0,0]';

numFlyerVertices = 24; % TODO derive from bodyVerticies

%% animation
for i=1:length(tAnim)

    % calculate angles for the orientation of the thrust vector before
    % applying quaternion rotation in the updated axis of rotation 
    % prop always rotates about center 
    % thrust vector along the axis of rotation, normal to the propeller
    % surface 
    
    %use unit vector as the direction for the thrust calculations  %use new alpha in quaternion to rotate prop
    %use rotation matrix to set the angle for axis of rotation used for quaternion below
    %angle for rotation is determined by the pitching induced from
    %eccentric thrust 
    axisRotateQuat= [(cos(zetaAnim3(i)/2)), sin(zetaAnim3(i)/2)*[0,0,1]];
    % Normalize Rotation Vector 
    axisRotateQuat = axisRotateQuat./sqrt((axisRotateQuat(1))^2+(axisRotateQuat(2))^2+(axisRotateQuat(3))^2);
    a1 = quaternion(axisRotateQuat);
    
    pitchRotate = [(cos(zetaAnim4(i)/2)), sin(zetaAnim4(i)/2)*[0,1,0]];
    pitchRotate = pitchRotate./sqrt((pitchRotate(1))^2+(pitchRotate(2))^2+(pitchRotate(3))^2);
    p1 = quaternion(pitchRotate);
    
    
    quatProp2World = [(cos(zetaAnim(i)/2)), sin(zetaAnim(i)/2)*[0 0 1]];
    
    % normalize quaternion 
    quatProp2World = quatProp2World./(sqrt(quatProp2World(1)^2+quatProp2World(2)^2+quatProp2World(3)^2+quatProp2World(4)^2));
    q1 = quaternion(quatProp2World);
    %convert quaternion to rotation matrix
    % TODO: Replace with inbuilt function quat2rotm
    rotm = quat2rotm(q1);
    rotm2 = quat2rotm(a1);
    rotm3 = quat2rotm(p1);
    
    % Homogenous transformation 
    propHT = [rotm, posWorld2Body; 0, 0, 0, 1];
    
    bodyHT = [rotm2, posWorld2Body; 0,0,0,1];
    flyerHT = [rotm3, [0;0;0]; 0, 0, 0, 1];
    for j = 1:numFlyerVertices 
        bodyVerticesNew(j,:) = (bodyHT*flyerHT*[bodyPatch.Vertices(j,:),1]');
        propVerticesNew(j,:) = propHT*flyerHT*[propPatch.Vertices(j,:),1]';
    end

    propVertices = propVerticesNew(:,1:3);
    bodyVertices = bodyVerticesNew(:,1:3);
    
    posWorld2Body = [zetaAnim5(i), zetaAnim6(i), zetaAnim2(i)]';
    %apply quaternion
    %vertex = rotatepoint(q1, vertexOrig);
    %draw figure
    figure(1)
    view(3) 
    propPatch.Vertices =  propVertices;
    bodyPatch.Vertices = bodyVertices;
   
    %pause length controls timestep
    pause(dt)
end
end