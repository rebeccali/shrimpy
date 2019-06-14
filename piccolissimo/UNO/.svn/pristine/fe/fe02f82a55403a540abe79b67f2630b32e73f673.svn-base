global beta_d chord_d aoa_d_out span_p span_d B_d Xbase nu_out
%% User changable params
plot_all = false;
modifications = {1,1,1,0,0,0};
reset_modifications = true;
modifications_orig = modifications;
test_states = [nan; .001; .001; nan; .001; .001; nan; nan; nan; .001; .001; nan; nan; nan; nan; nan; nan; nan];
% SetupFunctions = {@UnoOffset_Base_gemfan_10x4_5}; %{@UnoV3_gemfan_LL};%{@UnoV3_14Deg}; %{@UnoV3_gemfan_12x4_5}; %{@UnoV3_14Deg};
% SetupFunctions = {@UnoOffset_Base_gemfan_10x4_5, @UnoOffset_Base_gemfan_11x4_7, @UnoOffset_Base_gemfan_12x4_5};
SetupFunctions = {@UnoOffset_actual_gemfan_10x4_5, @UnoOffset_actual_gemfan_11x4_7, @UnoOffset_actual_gemfan_12x4_5};


% AOA
target_aoa_array = deg2rad([-2]); %deg2rad([-1.25 -1.5, -1.75, -2, -2.25, -2.5]);%deg2rad([-1.5 -2.25 -3 -3.75 -4.5]); %deg % fill color
rad_per_aoa = .9; % 1 means 1 rad off = 1 rad change
rad_per_span_array = 0;%deg2rad([0, -1, -2]); % marker
rad_per_nu = 0;%-.01; % 1 means 1m/s from average = 1 rad change
% Chord
target_rps_array = [60 20]; %line color range
m_per_rps = -.002;
m_per_nu = 0;
m_per_span_array = [.06];
chord_d_bulk_array = [.06]; % bulk chord change from default
% Params
iterations = 1;
% Limits
% chord_d_max = 1.25; % this value * the nominal chord_d
disc_solidity_max = .66; % Maximum allowable disc solidity, evenly distributed through span
chord_d_min = .03; % this value straight across
material_thickness = .04;

%% Init figures
figure(40);
clf;
hold all;
figure(41);
clf;
hold all;
figure(42);
clf;
hold all;
figure(43);
clf;
hold all;
figure(44);
clf;
hold all;
figure(45);
clf;
hold all;
figure(46);
clf;
hold all;
figure(47);
clf;
hold all;
figure(48);
clf;
hold all;
figure(49);
clf;
hold all;

%% Init variables
state_transitions = {};
target_aoa_out = {nan};
rad_per_span_out = {nan};
body_rate_out = {};
legend_out = {};
plot_legend_out = {};
m_per_span_out = {nan};
chord_d_bulk_out = {nan};

idx = 1;
maximin_eig = -inf;

aoa_color = [0 0 0];
chord_color = [0 0 0];
marker = 'x';
marker2 = 'x';

total_trials = length(rad_per_span_array)*length(target_aoa_array)*length(m_per_span_array)*length(chord_d_bulk_array)*iterations*length(SetupFunctions)

markers = {'o','s','d','p','h'};
speed_colors = varycolor(round(max(target_rps_array)-min(target_rps_array)));
aoa_colors = varycolor(length(target_aoa_array));
chord_colors = varycolor(length(chord_d_bulk_array));

%% Generate plot descriptions
rad_per_span_desc = {'Degrees per span marker:'; 'x = orig'};
for i = 1:length(rad_per_span_array)
    rad_per_span_desc{i+1,1} = [markers{i} ' = ' num2str(rad2deg(rad_per_span_array(i)))];
end

m_per_span_desc = {'Meters per span marker:'; 'x = orig'};
for i = 1:length(m_per_span_array)
    m_per_span_desc{i+1,1} = [markers{i} ' = ' num2str(m_per_span_array(i))];
end

speed_desc = {'Body speed, edge color RGBK:'; ['\color[rgb]{' num2str(speed_colors(1,:),'%f,%f,%f') '}' num2str(max(target_rps_array)) ' (rad/s)']; ['\color[rgb]{' num2str(speed_colors(end,:),'%f,%f,%f') '}' num2str(min(target_rps_array)) ' (rad/s)']};

aoa_desc = {'Average AOA, body color RGBK:'};
for i = 1:length(target_aoa_array)
    aoa_desc{i+1,1} = ['\color[rgb]{' num2str(aoa_colors(i,:),'%f,%f,%f') '}' num2str(rad2deg(target_aoa_array(i))) ' (deg)'];
end

chord_desc = {'Chord change, body color RGBK:'};
for i = 1:length(chord_d_bulk_array)
    chord_desc{i+1,1} = ['\color[rgb]{' num2str(chord_colors(i,:),'%f,%f,%f') '}' num2str(chord_d_bulk_array(i)) ' (m)'];
end

%% Optimization loops
for setup_function_idx = 1:length(SetupFunctions)
    setup_function = SetupFunctions(setup_function_idx);
    for rad_per_span_idx = 1:length(rad_per_span_array)
        rad_per_span = rad_per_span_array(rad_per_span_idx);
        for target_aoa_idx = 1:length(target_aoa_array)
            target_aoa = target_aoa_array(target_aoa_idx);
            for m_per_span_idx = 1:length(m_per_span_array)
                m_per_span = m_per_span_array(m_per_span_idx);
                for chord_d_bulk_idx = 1:length(chord_d_bulk_array)
                    chord_d_bulk = chord_d_bulk_array(chord_d_bulk_idx);
                    for itr = 1:iterations                 
                        %% Run trim the sim + find eigen values
                        if reset_modifications
                            modifications = modifications_orig;
                        end
                        state_transitions{end+1} = FlyerParamSweep( setup_function, modifications, test_states, idx-1 );
                        speed_color_idx = min(round(max(target_rps_array)-min(target_rps_array)),max(1,round(Xbase(9)-min(target_rps_array))));
                        speed_color = speed_colors(size(speed_colors,1) + 1 - speed_color_idx,:);
                        if idx == 1
    %                         chord_d_max = abs(2*pi*span_d./cos(beta_d)/B_d*disc_solidity_max);
                            chord_d_min = chord_d_min*ones(size(chord_d));
                            chord_d_orig = chord_d;
                            beta_d_mod = zeros(size(beta_d));
                            beta_d_orig = beta_d;
                        end

                        %% Print resulting body pitch angles
                        if itr == 1 || plot_all
                            plot_legend_out{end+1} = num2str(idx);

                            PlotOptimizations(idx, span_d, chord_d, beta_d, aoa_d_out, nu_out, speed_color, aoa_color, chord_color, marker, marker2, Xbase, state_transitions{idx})
                        end

                        %% Store results from this trial
                        target_aoa_out{end+1} = target_aoa;
                        rad_per_span_out{end+1} = rad_per_span;
    %                     target_rps_out{end+1} = target_rps;
                        body_rate_out{end+1} = Xbase(9);
                        legend_out{end+1} = num2str(idx);
                        m_per_span_out{end+1} = m_per_span;
                        chord_d_bulk_out{end+1} = chord_d_bulk;

                        aoa_color = aoa_colors(target_aoa_idx,:);
                        chord_color = chord_colors(chord_d_bulk_idx,:);
                        marker = markers{rad_per_span_idx};
                        marker2 = markers{m_per_span_idx};

                        %% Make modifications to chord and beta for next trial
                        mean_nu = sum((nu_out(2:end,end)' + nu_out(1:end-1,end)')/2 ... % average inflow between segments
                            .*(span_d(2:end).*span_d(2:end) - span_d(1:end-1).*span_d(1:end-1))*2*pi ... % area of a segment
                            /(span_d(end)*span_d(end)*2*pi)); % total area
                        if ~isempty(aoa_d_out)
                            beta_d_mod(chord_d ~= 0) = beta_d_mod(chord_d ~= 0) ... % last modification
                                + rad_per_aoa*(target_aoa - aoa_d_out(chord_d ~= 0,1,end)') ... % angle from target aoa
                                + rad_per_nu*(mean_nu-nu_out(chord_d ~= 0,end)') ... % error from uniform inflow
                                + rad_per_span*(span_d(chord_d ~= 0)/span_d(end)-.5);  % Linear lift distribution
                        else
                            beta_d_mod = zeros(size(beta_d_mod));
                        end

                        chord_d_mod(chord_d ~= 0) = chord_d_bulk + m_per_span*(span_d(chord_d ~= 0)/span_d(end)-.5);

                        chord_d_max_solidity = abs(2*span_d*sin((2*pi)./B_d*disc_solidity_max/2)./cos(beta_d_orig+beta_d_mod));
                        chord_d_max_thickness = material_thickness./abs(sin(beta_d_orig + beta_d_mod));
                        chord_d_max = min(chord_d_max_solidity,chord_d_max_thickness);

                        chord_d_mod(chord_d ~= 0) = min(chord_d_mod(chord_d ~= 0) ...
                                            , chord_d_max(chord_d ~= 0) - chord_d_orig(chord_d ~= 0));

                        chord_d_mod(chord_d ~= 0) = max(chord_d_mod(chord_d ~= 0) ...
                                            , chord_d_min(chord_d ~= 0) - chord_d_orig(chord_d ~= 0));

                        modifications(4) = {chord_d_mod};
                        modifications(5) = {beta_d_mod};

                        disp([num2str(idx/(total_trials+1)*100) '% complete']);
                        idx = idx+1;
                    end
                end
            end
        end
    end
end
%% Run trim the sim + find eigen values one last time
if reset_modifications
    modifications = modifications_orig;
end
state_transitions{end+1} = FlyerParamSweep( setup_function, modifications, test_states, idx-1 );
% target_aoa_out{end+1} = target_aoa;
% rad_per_span_out{end+1} = rad_per_span;
% target_rps_out{end+1} = target_rps;
body_rate_out{end+1} = Xbase(9);
legend_out{end+1} = num2str(idx);
plot_legend_out{end+1} = num2str(idx);
speed_color_idx = min(round(max(target_rps_array)-min(target_rps_array)),max(1,round(Xbase(9)-min(target_rps_array))));
speed_color = speed_colors(size(speed_colors,1) + 1 - speed_color_idx,:);
% aoa_color = aoa_colors(target_aoa_idx,:); % From previous run
% chord_color = chord_colors(chord_d_bulk_idx,:); % From previous run
PlotOptimizations(idx, span_d, chord_d, beta_d, aoa_d_out, nu_out, speed_color, aoa_color, chord_color, marker, marker2, Xbase, state_transitions{idx});

figure(40);
title('Stabilizer 1 aoa vs span (t(end))');
xlabel('Span (m)');
ylabel('AOA (deg)');
grid on;
% ylim([-20 20]);
xlim([0 max(span_p(end),span_d(end))]);
legend(plot_legend_out,'Interpreter','none','Location','west');
% axes('Position',[0 0 1 1],'Visible','off');
% text(0,.5,{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}});
annotation('textbox',[.25 .25 .3 .3],'String',{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}},'FitBoxToText','on');

figure(41);
title('Stabilizer beta vs span (before/after)')
xlabel('Span (m)');
ylabel('Beta (deg)');
grid on;
xlim([0 max(span_p(end),span_d(end))]);
legend(plot_legend_out,'Interpreter','none','Location','west');
% axes('Position',[0 0 1 1],'Visible','off');
% text(0,.5,{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}});
annotation('textbox',[.25 .25 .3 .3],'String',{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}},'FitBoxToText','on');

figure(42);
title('Propeller Speed')
xlabel('Iteration ()');
ylabel('Speed (rad/s)');
grid on;
legend(plot_legend_out,'Interpreter','none','Location','west');
% axes('Position',[0 0 1 1],'Visible','off');
% text(0,.5,{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}});
annotation('textbox',[.25 .25 .3 .3],'String',{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}},'FitBoxToText','on');

figure(43);
title('Body Speed')
xlabel('Iteration ()');
ylabel('Speed (rad/s)');
grid on;
legend(plot_legend_out,'Interpreter','none','Location','west');
% axes('Position',[0 0 1 1],'Visible','off');
% text(0,.5,{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}});
annotation('textbox',[.25 .25 .3 .3],'String',{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}},'FitBoxToText','on');

figure(44);
title('Inflow vs span (t(end))');
xlabel('Span (m)');
ylabel('Nu (m/s)');
grid on;
legend(plot_legend_out,'Interpreter','none','Location','west');
% axes('Position',[0 0 1 1],'Visible','off');
% text(0,.5,{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}});
annotation('textbox',[.25 .25 .3 .3],'String',{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}},'FitBoxToText','on');

figure(45);
title('Chord vs span');
xlabel('Span (m)');
ylabel('Chord (m)');
grid on;
legend(plot_legend_out,'Interpreter','none','Location','west');
% plot([span_d(end:-1:1) span_d(1:end)], [chord_d_min(end:-1:1) chord_d_max(1:end)]);
% legend({plot_legend_out{:}, 'limits'},'Interpreter','none','Location','west');
% axes('Position',[0 0 1 1],'Visible','off');
% text(0,.5,{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}});
annotation('textbox',[.25 .25 .3 .3],'String',{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}},'FitBoxToText','on');

state_transition = state_transitions{end};
temp_eig = complex(eig(state_transition{end}));

figure(46);
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
    temp_4_eig = sort(abs(xdata{idx}+ydata{idx}*1i));
    temp_4_eig = temp_4_eig(4);
    if temp_4_eig > maximin_eig
        maximin_eig = temp_4_eig;
    end
    
    temp_min_real_eig = sort(abs(xdata{idx}));
    min_real_eigs(length(xdata)-idx+1) = temp_min_real_eig(1);
end
[min_real_eigs_sorted, eigs_sorted_idx] = sort(min_real_eigs);
for i = 1:length(xdata)
    disp([plot_legend_out(eigs_sorted_idx(i)) num2str(min_real_eigs_sorted(i))]);
end
grid on;
axis square;
xlim([min(-maximin_eig,-.0001), 0]);
ylim([0, max(maximin_eig,.0001)]);
title('Eigen Values Zoomed');
xlabel('Real');
ylabel('Imaginary');
legend(plot_legend_out,'Interpreter','none','Location','west');
% axes('Position',[0 0 1 1],'Visible','off');
% text(0,.5,{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}});
annotation('textbox',[.25 .25 .3 .3],'String',{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}},'FitBoxToText','on');

figure(47);
grid on;
axis square;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('Eigen Values');
xlabel('Real');
ylabel('Imaginary');
legend(plot_legend_out,'Interpreter','none','Location','west');
% axes('Position',[0 0 1 1],'Visible','off');
% text(0,.5,{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}});
annotation('textbox',[.25 .25 .3 .3],'String',{rad_per_span_desc{:,1},speed_desc{:,1},aoa_desc{:,1}},'FitBoxToText','on');

figure(48);
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
    temp_4_eig = sort(abs(xdata{idx}+ydata{idx}*1i));
    temp_4_eig = temp_4_eig(4);
    if temp_4_eig > maximin_eig
        maximin_eig = temp_4_eig;
    end
end
grid on;
axis square;
xlim([min(-maximin_eig,-.0001), 0]);
ylim([0, max(maximin_eig,.0001)]);
title('Eigen Values Zoomed');
xlabel('Real');
ylabel('Imaginary');
legend(plot_legend_out,'Interpreter','none','Location','west');
% axes('Position',[0 0 1 1],'Visible','off');
% text(0,.5,{m_per_span_desc{:,1},speed_desc{:,1},chord_desc{:,1}});
annotation('textbox',[.25 .25 .3 .3],'String',{m_per_span_desc{:,1},speed_desc{:,1},chord_desc{:,1}},'FitBoxToText','on');

figure(49);
grid on;
axis square;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('Eigen Values');
xlabel('Real');
ylabel('Imaginary');
legend(plot_legend_out,'Interpreter','none','Location','west');
annotation('textbox',[.25 .25 .3 .3],'String',{m_per_span_desc{:,1},speed_desc{:,1},chord_desc{:,1}},'FitBoxToText','on');
% axes('Position',[0 0 1 1],'Visible','off');
% text(0,.5,{m_per_span_desc{:,1},speed_desc{:,1},chord_desc{:,1}});

beep;