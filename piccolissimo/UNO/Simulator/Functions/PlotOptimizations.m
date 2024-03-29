function PlotOptimizations(idx, span_d, chord_d, beta_d, aoa_d_out, nu_out, speed_color, aoa_color, chord_color, marker, marker2, Xbase, state_transition)
    figure(40);
    plot(span_d(chord_d ~= 0),rad2deg(aoa_d_out(chord_d ~= 0,1,end)) ...
        ,[marker '-'], 'LineWidth',2, ...
        'MarkerEdgeColor',speed_color, ...
        'MarkerFaceColor',aoa_color,'MarkerSize',10);
    hold all;

    figure(41);
    plot(span_d(chord_d ~= 0),rad2deg(beta_d(chord_d ~= 0)) ...
        ,[marker '-'], 'LineWidth',2, ...
        'MarkerEdgeColor',speed_color, ...
        'MarkerFaceColor',aoa_color,'MarkerSize',10);
    hold all;

    figure(42);
    plot(idx, Xbase(8)+Xbase(9) ...
        ,marker, 'LineWidth',2, ...
        'MarkerEdgeColor',speed_color, ...
        'MarkerFaceColor',aoa_color,'MarkerSize',10);
    hold all;

    figure(43);
    plot(idx, Xbase(9) ...
        ,marker, 'LineWidth',2, ...
        'MarkerEdgeColor',speed_color, ...
        'MarkerFaceColor',aoa_color,'MarkerSize',10);
    hold all;

    figure(44);
    plot(span_d,nu_out(:,end) ...
        ,[marker '-'], 'LineWidth',2, ...
        'MarkerEdgeColor',speed_color, ...
        'MarkerFaceColor',aoa_color,'MarkerSize',10);
    hold all;

    figure(45);
    plot(span_d,chord_d ...
        ,[marker '-'], 'LineWidth',2, ...
        'MarkerEdgeColor',speed_color, ...
        'MarkerFaceColor',aoa_color,'MarkerSize',10);
    hold all;

    temp_eig = complex(eig(state_transition{end}));
    figure(46);
    plot(temp_eig ...
        ,marker, 'LineWidth',2, ...
        'MarkerEdgeColor',speed_color, ...
        'MarkerFaceColor',aoa_color,'MarkerSize',10);
    hold all;

    figure(47);
    plot(temp_eig ...
        ,marker, 'LineWidth',2, ...
        'MarkerEdgeColor',speed_color, ...
        'MarkerFaceColor',aoa_color,'MarkerSize',10);
    hold all;

    figure(48);
    plot(temp_eig ...
        ,marker2, 'LineWidth',2, ...
        'MarkerEdgeColor',speed_color, ...
        'MarkerFaceColor',chord_color,'MarkerSize',10);
    hold all;

    figure(49);
    plot(temp_eig ...
        ,marker2, 'LineWidth',2, ...
        'MarkerEdgeColor',speed_color, ...
        'MarkerFaceColor',chord_color,'MarkerSize',10);
    hold all;
end