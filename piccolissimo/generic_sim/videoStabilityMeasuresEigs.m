function [  ] = videoStabilityMeasuresEigs(testVMag, testOmgMag)
global h_d Ir_r Is_s rotorErrorAngle
xlimit = [-1.5, 1.5];
ylimit = [-100 100];
h_d_base = h_d;
h_ds = h_d-.015:.002:h_d+.015;
Ir_r_base = Ir_r;
Is_s_base = Is_s;
Ir_r_ratios = .5:.05:1;
Is_s_ratios = .5:.05:1;
% rotorErrorAngle = eye(3);
angle_values = 0:10;
writerObj = VideoWriter('stabilityMeasuresEigs.mp4','MPEG-4');
writerObj.FrameRate = 15;
open(writerObj);
figure(20);
close gcf;
figure(20);
% hold on;
grid on;
for i = 1:length(h_ds)
    h_d = h_ds(i);
    [eigVals rh accels accels2 forces] = StabilityMeasures(testVMag,testOmgMag);
    figure(20);
    plot(real(eigVals),imag(eigVals),'*','Color',[i/length(h_ds)/2+.5, 0, 0]) %   Plot real and imaginary parts
    xlim(xlimit);
    ylim(ylimit);
    grid on;
    set(gca,'xtick',[xlimit(1) 0 xlimit(2)]);
    set(gca,'ytick',[ylimit(1) 0 ylimit(2)]);
    set(gcf,'Color',[1 1 1]);
    title(['h_d = ' num2str(h_d*1000) 'mm']);
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end
h_d = h_d_base;

% for i = 1:length(Ir_r_ratios)
%     Ir_r = Ir_r_base.*Ir_r_ratios(i);
%     [eigVals rh accels accels2 forces] = StabilityMeasures(testVMag,testOmgMag);
%     figure(20);
%     plot(real(eigVals),imag(eigVals),'*','Color',[0, i/(2*length(Ir_r_ratios))/2+.5, 0]) %   Plot real and imaginary parts
%     xlim(xlimit);
%     ylim(ylimit);
%     grid on;
%     title('Increasing I_R');
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
% end
% for i = length(Ir_r_ratios):-1:1
%     Ir_r = Ir_r_base./Ir_r_ratios(i);
%     [eigVals rh accels accels2 forces] = StabilityMeasures(testVMag,testOmgMag);
%     figure(20);
%     plot(real(eigVals),imag(eigVals),'*','Color',[0, (length(Ir_r_ratios)-i+1)/(2*length(Ir_r_ratios))/2+.75, 0]) %   Plot real and imaginary parts
%     xlim(xlimit);
%     ylim(ylimit);
%     grid on;
%     title('Increasing I_R');
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
% end
% Ir_r = Ir_r_base;
% 
% for i = 1:length(Is_s_ratios)
%     Is_s = Is_s_base.*Is_s_ratios(i);
%     [eigVals rh accels accels2 forces] = StabilityMeasures(testVMag,testOmgMag);
%     figure(20);
%     plot(real(eigVals),imag(eigVals),'*','Color',[0, 0, i/(2*length(Is_s_ratios))/2+.5]) %   Plot real and imaginary parts
%     xlim(xlimit);
%     ylim(ylimit);
%     grid on;
%     title('Increasing I_S');
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
% end
% for i = length(Is_s_ratios):-1:1
%     Is_s = Is_s_base./Is_s_ratios(i);
%     [eigVals rh accels accels2 forces] = StabilityMeasures(testVMag,testOmgMag);
%     figure(20);
%     plot(real(eigVals),imag(eigVals),'*','Color',[0, 0, (length(Is_s_ratios)-i+1)/(2*length(Is_s_ratios))/2+.75]) %   Plot real and imaginary parts
%     xlim(xlimit);
%     ylim(ylimit);
%     grid on;
%     title('Increasing I_S');
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
% end
% Is_s = Is_s_base;
% 
% for i = 1:length(angle_values)
%     rotorErrorxAngle = deg2rad(angle_values(i));
%     rotorErroryAngle = deg2rad(0);
%     rotorErrorAngle = [1,0,0;0,cos(rotorErrorxAngle),-sin(rotorErrorxAngle);0,sin(rotorErrorxAngle),cos(rotorErrorxAngle)]*[cos(rotorErroryAngle),0,sin(rotorErroryAngle);0,1,0;-sin(rotorErroryAngle),0,cos(rotorErroryAngle)];
%     [eigVals rh accels accels2 forces] = StabilityMeasures(testVMag,testOmgMag);
%     figure(20);
%     plot(real(eigVals),imag(eigVals),'*','Color',[i/length(angle_values)/2+.5, i/length(angle_values)/2+.5, 0]) %   Plot real and imaginary parts
%     xlim(xlimit);
%     ylim(ylimit);
%     grid on;
%     title('Increasing error');
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
% end
% rotorErrorAngle = eye(3);

close(writerObj);
end
