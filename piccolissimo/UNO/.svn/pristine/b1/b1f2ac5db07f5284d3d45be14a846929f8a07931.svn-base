function [ state_transitions ] = FlyerParamSweep( SetupFunction, modifications, test_states, add_param )
%   [Is_s, R_d, B_d, chord_d, beta_d, beta_p]
    global Xbase
    if nargin <= 3 || ~add_param
        figure(30);
        clf;
    end
    maximin_eig = -inf;
    state_transition_filter = ones(6);
%     state_transition_filter = [1 0 0 0 0 1; ...
%                                0 1 0 0 1 0; ...
%                                1 0 1 1 0 0; ...
%                                0 1 1 1 0 0; ...
%                                0 0 1 0 0 0; ...
%                                0 0 0 1 0 0];
    
    for idx = 1:length(SetupFunction)
        for k=1:size(modifications,1)
%             keyboard
%             disp([modifications{k,:}]);
            SetupFunction{idx}();
            current_ind = (idx-1)*size(modifications,1) + k;
            ModifySetup(modifications(k,:));
            [state_transitions{current_ind}, ~] = GenerateStateTransition(test_states);
            %% Clean up state transition to only have a, b, d, e, g
            state_transitions{current_ind} = state_transitions{current_ind}.*state_transition_filter;
            if nargin >3
                state_labels{current_ind} = num2str(add_param+1);
            else
                state_labels{current_ind} = [func2str(SetupFunction{idx}) num2str([modifications{k,:}],'% +6.2f')];
            end
            disp(Xbase');
%             damp(state_transitions{current_ind})
%             [V,D] = eig(state_transitions{current_ind})
            figure(30); 
            temp_eig = eig(state_transitions{current_ind});
            
            subplot(1,2,1);
            plot(complex(temp_eig),'x','LineWidth',2,'MarkerSize',10,'DisplayName',state_labels{current_ind});
            hold all;
            legend('-DynamicLegend');
            
%             figure(16);
            subplot(1,2,2);
            plot(complex(temp_eig),'x','LineWidth',2,'MarkerSize',10,'DisplayName',state_labels{current_ind});
            hold all;
            legend('-DynamicLegend');
        end
    end
    
    figure(30);
    subplot(1,2,1);
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    dataObjs = get(ax, 'Children');
    xdata = get(dataObjs, 'XData');
    ydata = get(dataObjs, 'YData');
    if isa(xdata,'double')
        xdata = {xdata};
        ydata = {ydata};
    end
    for idx = 1:length(xdata)
        temp_min_eig = sort(abs(xdata{idx}+ydata{idx}*1i));
        temp_min_eig = temp_min_eig(4);
        if temp_min_eig > maximin_eig
            maximin_eig = temp_min_eig;
        end
    end
    grid on;
    axis square;
    xlim([min(-maximin_eig,-.0001), 0]);
    ylim([0, max(maximin_eig,.0001)]);
    
    subplot(1,2,2);
    grid on;
    axis square;
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    
%     beep
end

function ModifySetup(modifications)
    global Is_s R_d1 R_d2 B_d beta_d beta_p chord_d

    Is_s = modifications{1}*Is_s;
    R_d1 = modifications{2}*R_d1;
    R_d2 = R_d1;
    B_d = round(modifications{3}*B_d);
    chord_d = modifications{4}+ chord_d;
    beta_d = modifications{5} + beta_d;
    beta_p = modifications{6} + beta_p;
end