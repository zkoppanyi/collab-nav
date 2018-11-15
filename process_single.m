clear variables;

%% Settings 
simulation_settings;

%gps_error_tests = [0.1 0.5 1.0 5.0 10.0];
    
load('problem');
clear agents


% Clear persistent variables
clear SelfishAgent

%% Initialization
roads = problem.roads;
% plot(problem.vehicles(:, 2), problem.vehicles(:, 3), 'k.')
% axis equal;

% Get timestamps
ts = problem.vehicles(:,1);
ts = unique(ts);

%form_t = 0; to_t = 1000;
form_t = 100; to_t = 110;

dxy_all = [];

% Discover how many vehicles we have in the give time interval
veh_ids_total = [];
for t = form_t : dt*5 : to_t
    idx = find( abs(problem.vehicles(:, 1) - t) < dt/2 );
    epoch = problem.vehicles(idx, :);    
    veh_ids = unique(epoch(:,2));
    veh_ids_total = [veh_ids_total; veh_ids]; 
end
veh_ids_total = unique(veh_ids_total);
veh_ids_total = [10 16 8 20 13 6 17];  % fixed ids
%agents{length(veh_ids_total)} = {};
agents{250} = {};

figure(1);
set(gcf,'Position',[50 200 500 500])
figure(2); clf; hold on;                
set(gcf,'Position',[600 200 800 700])

%% Algorithm

% Init filter
n_veh = length(veh_ids_total);
n_states = 4;
x = zeros(n_veh*n_states, 1);
P = eye(n_veh*n_states, n_veh*n_states);

for t = form_t : dt : to_t

    fprintf("t= %.1f\n", t); 
    
    idx = find( abs(problem.vehicles(:, 1) - t) < dt/2 );
    epoch = problem.vehicles(idx, :);
    
    if isempty(epoch)
        continue;
    end
    
    % Find trajectory within the AOI
    epoch = epoch(and(and(and(epoch(:,3) > system_setting.AOI (1,1), epoch(:,3) < system_setting.AOI (1,2)), epoch(:,4) > system_setting.AOI (2,1)), epoch(:,4) < system_setting.AOI (2,2)), :);   
           
    % Vehicles at the current epoch
    veh_ids = unique(epoch(:,2));       
%     if isempty(find(veh_ids == sel_veh_id, 1))
%         continue;
%     end

    %Fixed number of vehicles
    rm_idx = [];
    for i = 1 : length(veh_ids)
        if isempty(find(veh_ids(i) == veh_ids_total))
            rm_idx = [rm_idx; i];
        end
    end
    veh_ids(rm_idx) = [];
    
     if isempty(veh_ids)
         continue;
     end

    %% Internal observations
    
    for i = 1 : length(veh_ids)
        
        veh_id = veh_ids(i);      
                
        % Get internal observations
        inter_obs = epoch(epoch(:,2) == veh_id, :);
                 
        % Get internal observation
        xy = inter_obs(3:4)';
        b  = inter_obs(5)/180*pi;
        v  = inter_obs(6);
        upt = [xy', v, b, t, 0];
        
        % Initialization
        if isempty(agents{veh_id})
           x_init = [xy; v; b];
           P_init = eye(4,4);  
           agents{veh_id} = SelfishAgent(veh_id, x_init, P_init, t, n_veh, system_setting);
        end        
            
        % Build prediction: F, f, Q
        agents{veh_id}.build_predict();     
     
        % Measurement\observation        
        if mod(t, 0.2) == 0
            agents{veh_id}.build_int_update('GPS+IMU', upt);
        else
            agents{veh_id}.build_int_update('IMU', upt);            
        end
        
    % Extended Kalman filter
    [x, P] = ekf_predict(agents{veh_id}.f, agents{veh_id}.F, agents{veh_id}.x, agents{veh_id}.P, agents{veh_id}.Q);
    [x, P] = ekf_update(x, agents{veh_id}.z_int, agents{veh_id}.h_int, agents{veh_id}.H_int, P, agents{veh_id}.R_int);
    
    % Unscented Kalman filter
%     [x, P, X1, X2] = ukf_predict(f, x, P, Q);
%     [x, P]  = ukf_update(x, z, h, P, R, X1, X2);

    agents{veh_id}.apply_update(x, P);

    end            
       
    %% Visualization and analytics
    
    % Vehicles positions and communication links
    figure(1); clf; hold on;
    title(sprintf('t= %.1f Acc: %.2f', t, system_setting.sigma_GPS));
    plot(roads(:, 1), roads(:, 2), 'k.');    
    for i = 1 : length(agents)        
        if isempty(agents{i}), continue; end
        
        x_loc     = agents{i}.x;
        xy_loc_gt = agents{i}.gt(end, :);
        z_loc     = agents{i}.z_int_hist(end, :);

        plot(x_loc(1), x_loc(2), 'b.', 'MarkerSize', 15);
        plot(xy_loc_gt(1), xy_loc_gt(2), 'g.', 'MarkerSize', 15);  
        plot(z_loc(1), z_loc(2), 'r.', 'MarkerSize', 15); 

        if agents{i}.id == sel_veh_id
            text(xy_loc_gt(1), xy_loc_gt(2)+20, sprintf('%i',  agents{i}.id), 'BackgroundColor', 'w', 'Color', 'r');
        else
            text(xy_loc_gt(1), xy_loc_gt(2)+20, sprintf('%i',  agents{i}.id), 'BackgroundColor', 'w');
        end        
        
    end

    grid on; axis equal;    
    xlim(system_setting.AOI(1,1:2)); ylim(system_setting.AOI(2,1:2));       
    set(gca, 'FontSize', 12);    
    xlabel('[m]'); ylabel('[m]');
    
    % Selected vehicle's internal states    
    agent = agents{sel_veh_id};
    if ~isempty(agent)
        n = size(agent.x_hist, 1);
        dxy = (agent.x_hist(end, 1:4)-agent.gt(end, 1:4))';

        figure(2); clf; hold on;                

        subplot(4, 1, 1); hold on; 
        plot(1:n, agent.x_hist(:, 1), 'b.-' );
        plot(1:n, agent.gt(:, 1), 'g.-' );
        title(sprintf('[%i] X= %.3f', sel_veh_id, dxy(1)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 2); hold on; 
        plot(1:n, agent.x_hist(:, 2), 'b.-' );
        plot(1:n, agent.gt(:, 2), 'g.-' );
        title(sprintf('Y= %.3f', dxy(2)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 3); hold on; 
        plot(1:n, agent.x_hist(:, 3), 'b.-' );
        plot(1:n, agent.gt(:, 3), 'g.-' );
        title(sprintf('v= %.3f', dxy(3)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 4); hold on; 
        plot(1:n, agent.x_hist(:, 4)/pi*180, 'b.-' );
        plot(1:n, agent.gt(:, 4)/pi*180, 'g.-' );
        title(sprintf('bearing= %.3f', dxy(4)/pi*180));  
        set(gca, 'FontSize', 10);       
    end
        
    
    % Show P just those elements that are available
    figure(3);
    show_sparisty_pattern(P, 3);
    
    % Compate trajectries to ground truth   
    fprintf('Time: %.1f\n', t);    
    for i = 1 : length(veh_ids)
        if ~isempty(agents{veh_ids(i)})
            dxy = agents{veh_ids(i)}.gt(:, 1:2) - agents{veh_ids(i)}.x_hist(:, 1:2);
            dxy = sqrt(dxy(:, 1).^2 + dxy(:, 2).^2);
            xy_rmse = sqrt(sum(dxy.^2) / length(dxy));
            xy_std = std(dxy);
            fprintf('#%i RMSE: [%.3f]  STD: [%.3f]\n', veh_ids(i), xy_rmse, xy_std);
            dxy_all  = [dxy_all; dxy(end)];
        end
    end
    fprintf('TOTAL MEAN: [%.3f] MEDIAN: [%.3f]\n', mean(dxy_all), median(dxy_all));
    
    pause(viz_frame_rate);
end

sol.agents = agents;
sol.system_setting = system_setting;
sol.problem = problem;
save(['results\' sprintf('solution_single_%i', simulation_scenario)], 'sol')


