function plotFlyer(args, tout, Xout, velocity_world, orientation)
global thrust_p_out pwm_out v_clamp r_b rho chord_d chord_p aoa_p_out aoa_d_out tau_p_out thrust_d_out RP_tau_p_out RP_tau_d_out i_m_out nu_out Ir_r Is_s span_d span_p I_tot K_t_m  v_out 


        figure(1);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,Xout(:,15),'color','r'); %x
        hold on;
        title('Position in World frame');
        plot(tout,Xout(:,16),'color','g'); %y
        plot(tout,Xout(:,17),'color','b'); %z

        figure(2); 
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,Xout(:,10),'color','r'); %phi
        hold on;
        title('Angle From World to Body');
        plot(tout,Xout(:,11),'color','g'); %theta
        plot(tout,Xout(:,12),'color','b'); %psi
    %     plot(tout,Xout(:,13),'color','m'); %rotor psi
    %     plot(tout,Xout(:,14),'color','c'); %stator psi

        figure(3);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,Xout(:,2),'color','r');%u
        hold on;
        title('Velocity in Flyer Frame');
        plot(tout,Xout(:,3),'color','g');%v
        plot(tout,Xout(:,4),'color','b');%w
        
        figure(33);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,velocity_world(:,1),'color','r');
        hold on;
        title('Velocity in World Frame');
        plot(tout,velocity_world(:,2),'color','g');
        plot(tout,velocity_world(:,3),'color','b');

        figure(4);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,Xout(:,5),'color','r');%p
        hold on;
        title('Angular Velocity in Flyer Frame');
        plot(tout,Xout(:,6),'color','g');%q
        plot(tout,Xout(:,7),'color','b');%r
        plot(tout,Xout(:,8),'color','c');%omg_r
        plot(tout,Xout(:,9),'color','y');%omg_b

        figure(5);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
    %     plot(tout,Xout(:,1)); %nu
        plot(tout(1:size(nu_out,2)),nu_out);
        hold all;
        title('Nu');
        plot(tout(1:size(nu_out,2)),mean(nu_out,1),'LineWidth',2);
        nu_legend = cellstr(num2str(span_d(:)));
        nu_legend{end+1} = 'Mean';
        legend(nu_legend);

        figure(6);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout,orientation(:,1),'color','r');%rotation
        hold on;
        title('Orientation of Spin Axis');
        plot(tout,orientation(:,2),'color','g');%rotation
        plot(tout,orientation(:,3),'color','b');%rotation

        figure(7);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
    %     plot(tout,Xout(:,18),'color','r');
        plot(tout(1:length(i_m_out)),i_m_out,'color','r');
        hold on;
        title('Current and Voltage');
        plot(tout(1:length(v_out)),v_out,'color','g');
        plot(tout(1:length(v_out)),v_clamp-pwm_out.*i_m_out*r_b,'color','b');
        legend('i','v_m','v_b');

%         figure(8);
%         hold off;
%         plot(tout(1:end-1),psi_des,'color','r');
%         hold on;
%         title('Commanded and error psi');
%         plot(tout,errpsis,'color','g');
%         mean(errpsis)
%         mean(errpsis(10000:end))
%         mean(errpsis(25000:end))
%         
% Fig 9 stuff here
        figure(11); % Angular momentum
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        z = zeros(length(Xout(:,8)),2);
        plot(tout,Ir_r*(Xout(:,5:7)+[z, Xout(:,8)])',tout,Is_s*(Xout(:,5:7)+[z, Xout(:,9)])',tout,I_tot*Xout(:,5:7)');
        legend('Rx','Ry','Rz','Sx','Sy','Sz','Bx','By','Bz');
        title('Momentum');

        figure(12); % Thrust sources
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:length(thrust_p_out)),thrust_p_out, tout(1:length(thrust_d_out)), thrust_d_out, tout(1:length(thrust_d_out)), thrust_d_out+thrust_p_out);
        legend('Propeller','Stabilizer','Net');
        title('Thrust');

        figure(13); % Torque sources
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:length(RP_tau_p_out)),RP_tau_p_out, tout(1:length(RP_tau_d_out)), RP_tau_d_out, tout(1:length(RP_tau_d_out)),RP_tau_p_out+RP_tau_d_out);
        legend('Propeller','Stabilizer','Net');
        title('Torque');

        figure(14); % Motor torque
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:length(i_m_out)),i_m_out*K_t_m);
        title('Motor Torque');

        figure(15); % Prop torque
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:length(i_m_out)),tau_p_out);
        title('Propeller torque directions');

        figure(16);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:size(aoa_p_out,3)),reshape(aoa_p_out(:,1,:),[size(aoa_p_out,1),size(aoa_p_out,3)]));
        title('Propeller 1 aoa vs time');

        figure(17);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(tout(1:size(aoa_d_out,3)),reshape(aoa_d_out(:,1,:),[size(aoa_d_out,1),size(aoa_d_out,3)]));
        title('Drag 1 aoa vs time');

        figure(18);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(span_p(chord_p ~= 0),rad2deg(aoa_p_out(chord_p ~= 0,1,end)),'.-')
        title('Propeller 1 aoa vs span (t(end))');
        xlabel('Span (m)');
        ylabel('AOA (deg)');
        grid on;
        ylim([-20 20]);
        xlim([0 max(span_p(end),span_d(end))]);

        figure(19);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(span_d(chord_d ~= 0),rad2deg(aoa_d_out(chord_d ~= 0,1,end)),'.-');
        title('Stabilizer 1 aoa vs span (t(end))');
        xlabel('Span (m)');
        ylabel('AOA (deg)');
        grid on;
        ylim([-20 20]);
        xlim([0 max(span_p(end),span_d(end))]);

        figure(20);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        plot(span_d,nu_out(:,end),'.-');
        title('Inflow vs span (t(end))');
        xlabel('Span (m)');
        ylabel('Nu (m/s)');
        grid on;

        figure(21);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        kinematic_viscosity = 1.846e-5/rho; % ~1.568e-5?
        plot(span_d(chord_d ~= 0), abs(Xout(end,9))*span_d(chord_d ~= 0).*chord_d(chord_d ~= 0)/kinematic_viscosity);
        title('Stabilizer Reynolds Number vs Span (t(end))')
        xlabel('Span (m)');
        ylabel('Reynolds Number ()');
        grid on;

        figure(22);
        ClearPlot(nnz(strcmpi(args,'plot_hold')) == 0);
        kinematic_viscosity = 1.846e-5/rho; % ~1.568e-5?
        plot(span_p(chord_p ~= 0), abs(Xout(end,8))*span_p(chord_p ~= 0).*chord_p(chord_p ~= 0)/kinematic_viscosity);
        title('Propeller Reynolds Number vs Span (t(end))')
        xlabel('Span (m)');
        ylabel('Reynolds Number ()');
        grid on;
        
        % Save all figures
        folderName = './generatedFigures';   % Your destination folder
        display(strcat('Saving all files to ', folderName));
        figList = findobj(allchild(0), 'flat', 'Type', 'figure');
        for iFig = 1:length(figList)
          figHandle = figList(iFig);
          figTitle = figHandle.Children(1).Title.String;
          if isempty(figTitle)
              figTitle = sprintf('%d_plot', iFig);
          end 
          figName   = strcat(figTitle, '.fig');
          savefig(figHandle, fullfile(folderName, figName));
        end
        close all;
        

end