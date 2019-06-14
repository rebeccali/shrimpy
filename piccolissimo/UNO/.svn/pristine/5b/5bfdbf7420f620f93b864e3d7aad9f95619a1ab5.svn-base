function [state_transition, inds_performed] = GenerateStateTransition(disturbances)
    global Xbase
    
    % Ensure we're trimmed and all globals initialized
    [ Xbase, found_trim ] = FlyerTrim(50);
    if found_trim
        Xbase([2,3,4,5,6,7,10,11,12,13,14,15,16,17]) = 0; % clear things that should be blank
        if size(disturbances) ~= size(Xbase)
            error('Xbase and disturbances sizes do not match');
        end

        inds_performed = find(~isnan(disturbances));

        for idx = 1:length(inds_performed)
            X_test = Xbase;
            X_test(inds_performed(idx)) = X_test(inds_performed(idx)) + disturbances(inds_performed(idx));
            dx_test = StateTransitionTest(X_test); % find dx
            state_transition(idx,:) = dx_test/disturbances(inds_performed(idx)); % Normalize by disturbance
        end
        % Trim state_transition to just inds that we care about
        state_transition = state_transition(:,inds_performed)';
        
    else
        inds_performed = zeros(size(disturbances));
        state_transition = zeros(6);
    end
end

function [means] = StateTransitionTest(X)
    %setup blade orientations
    rotation_positions = 24;
    rotation_step = pi/rotation_positions;
    k = 1;
    for i = 0:(rotation_positions-1)
        for j = 0:(rotation_positions-1)
            X_test = X;
            X_test(13) = X_test(13) + i*rotation_step;
            X_test(14) = X_test(14) + j*rotation_step;
            %run test
            [dX(k,:)] = FlyerODE(0,X_test);
            k = k+1;
        end
    end
    %sum tests
    means = mean(dX,1);
end