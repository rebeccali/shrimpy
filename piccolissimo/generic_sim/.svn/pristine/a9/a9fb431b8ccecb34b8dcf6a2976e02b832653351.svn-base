function [ state_transitions ] = FlyerParamSweep( SetupFunction, modifications, test_states )
%UNTITLED2 Summary of this function goes here
%   [Is_s, R_d, B_d, chord_d, beta_d, beta_p]
    global Xbase
    figure(15);
    clf;
    
    for i = 1:length(SetupFunction)
        for k=1:size(modifications,1)
            disp(modifications(k,:));
            SetupFunction{i}();
            ModifySetup(modifications(k,:));
            [state_transitions{k}, ~] = GenerateStateTransition(test_states);
            state_labels{(i-1)*size(modifications,1) + k} = [func2str(SetupFunction{i}) num2str(modifications(k,:),'% +6.2f')];
            disp(Xbase');
            damp(state_transitions{k})
            figure(15); 
            plot(eig(state_transitions{k}),'x','LineWidth',2,'MarkerSize',10);
            hold all;
        end
    end
    
    figure(15);
    legend(state_labels,'Interpreter','none','Location','east');
    grid on;
    axis square;
    xlim([-2.5 0]);
    ylim([0 2.5]);
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
end

function ModifySetup(modifications)
    global Is_s R_d1 R_d2 B_d beta_d beta_p chord_d

    Is_s = modifications(1)*Is_s;
    R_d1 = modifications(2)*R_d1;
    R_d2 = R_d1;
    B_d = round(modifications(3)*B_d);
    chord_d = modifications(4)*chord_d;
    beta_d = modifications(5) + beta_d;
    beta_p = modifications(6) + beta_p;
end