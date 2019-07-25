function plotFlyer(t, states)
%% Plot the flyer states
% TODO: add titles, standard color interfacing

figure(2);
plot(t,states(:,1),'b');
hold on;
plot(t,states(:,2),'r');
legend('prop angle (rad)', 'prop angular rate (rad/s)');

figure(3);
plot(t, states(:,5),'g');
hold on;
plot(t, states(:,6),'k');
legend('bdy angle (rad)', 'body angle rate (rad/s)');

figure(4)
plot(t, states(:,7),'b');
hold on;
plot(t, states(:,8),'k');
legend('pitch angle (rad)', 'pitch angle rate (rad/s)');

figure(5)
plot(t, states(:,9),'b');
hold on;
plot(t, states(:,10),'k');
legend('X position (m)', 'X velocity(m/s)');

figure(6)
plot(t, states(:,11),'b');
hold on;
plot(t, states(:,12),'k');
legend('Y position (m)', 'Y velocity(m/s)');
end