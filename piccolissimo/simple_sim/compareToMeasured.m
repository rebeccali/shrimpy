

voltage = [0
1.385
2.127
3.051
3.161];
current = [0
0.1512
0.2721
0.4392
0.4701];
speed = [0
23430
33090
42920
44510]*2*pi/60;
thrust = [0
0.01432
0.02922
0.05029
0.05405];

masses = .000:.001:.005;

for i = 1:length(masses)
    [ t,X, t_out, i_out, v_out, tau_out, thrust_out ] = simulateMotorPulse(masses(i));
    voltage_sim(i) = v_out(1);
    current_sim(i) = i_out(1);
    speed_sim(i) = X(1,5);
    thrust_sim(i) = thrust_out(1);
end
figure(1); clf;
plot(voltage_sim,current_sim,voltage,current);
title('Current vs Voltage');
legend('Sim','Act');
figure(2);clf;
plot(voltage_sim,speed_sim,voltage,speed);
title('Speed vs Voltage');
legend('Sim','Act');
figure(3);clf;
plot(voltage_sim,thrust_sim,voltage,thrust);
title('Thrust vs Voltage');
legend('Sim','Act');
figure(4);clf;
plot(speed_sim,current_sim,speed,current);
title('Current vs Speed');
legend('Sim','Act');
figure(5);clf;
plot(speed_sim,thrust_sim,speed,thrust);
title('Thrust vs Speed');
legend('Sim','Act');
figure(6);clf;
plot(current_sim,thrust_sim,current,thrust);
title('Thrust vs Current');
legend('Sim','Act');
