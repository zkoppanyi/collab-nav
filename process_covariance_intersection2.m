clear variables;

% Clear persistent variables
clear CoopAgent

load('problem');

%% Related materials:

% [1] Overview of federated filtering: "Federated Filtering Revisited: New Directions to Distributed Systems Estimation and Filtering – a Case Study"
% [2] Original F-EKF paper (1988): https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=195473
% [3] Investigation on information sharing factor: https://pubs.acs.org/doi/pdf/10.1021/ie0511175
% [4] Proposed information sharing factor based on local filters covariances: https://pdfs.semanticscholar.org/5eea/21a9db2448c07234b99c94ed6015c27d643c.pdf
% [5] Information-sharing factor based on median and predicated states: "Federated Filtering Revisited: New Directions to Distributed Systems Estimation and Filtering – a Case Study"
% [6] Original covariance intersection paper: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=609105



%% Settings 

system_setting.AOI = [750 1100; 1250 1600];                                 % Area of interest
viz_frame_rate = 0.1;                                                       % Frame rate for visualization
dt = 0.05;

% Error characterization
simulation_settings;

%% Initialization
roads = problem.roads;
% plot(problem.vehicles(:, 2), problem.vehicles(:, 3), 'k.')
% axis equal;
sel_veh_id  = 10; 

% Get timestamps
ts = problem.vehicles(:,1);
ts = unique(ts);
form_t = 100; to_t = 130;
dxy_all = [];

%problem.vehicles = problem.vehicles(problem.vehicles(:, 2) ~= 9, :);

% Discover how many vehicles we have in the give time interval
veh_ids_total = [];
for t = form_t : dt*5 : to_t
    idx = find( abs(problem.vehicles(:, 1) - t) < dt/2 );
    epoch = problem.vehicles(idx, :);    
    veh_ids = unique(epoch(:,2));
    veh_ids_total = [veh_ids_total; veh_ids]; 
end
veh_ids_total = unique(veh_ids_total);
%veh_ids_total = [10 16 8 20 13 6 17];
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
x_m = zeros(n_veh*n_states, 1);
P_m = eye(n_veh*n_states, n_veh*n_states);

init_rewind = form_t;
init_i = 0;
init_time = 0.5/dt;
is_init = 0;
weights = [];
t = form_t;
iter_i = 0;

while 1

    if t >= to_t
        break;
    end
    
    if init_i < init_time
        t = init_rewind;
        init_i = init_i + 1;
        is_init = 1;
    else
        is_init = 0;
        t = t + dt;
        iter_i = iter_i + 1;
    end
    
    fprintf("Time: %.3f\n", t);      
    
    idx = find( abs(problem.vehicles(:, 1) - t) < dt/2 );
    epoch = problem.vehicles(idx, :);
    
    if isempty(epoch)
        continue;
    end
    
    % Find trajectory within the AOI
    epoch = epoch(and(and(and(epoch(:,3) > system_setting.AOI (1,1), epoch(:,3) < system_setting.AOI (1,2)), epoch(:,4) > system_setting.AOI (2,1)), epoch(:,4) < system_setting.AOI (2,2)), :);   
           
    % Vehicles at the current epoch
    veh_ids = unique(epoch(:,2)); 
    
    % Fixed number of vehicles
    rm_idx = [];
    for i = 1 : length(veh_ids)
        if isempty(find(veh_ids(i) == veh_ids_total))
            rm_idx = [rm_idx; i];
        end
    end
    veh_ids(rm_idx) = [];
    
%     if isempty(find(veh_ids == sel_veh_id, 1))
%         continue;
%     end
    
    % New cycle; reset internal variables
     for i = 1 : length(agents)
        if isempty(agents{i}), continue; end            
        agents{i}.reset();
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
        upt = [xy', v, b, t, is_init];  
        
        % Initialization
        if isempty(agents{veh_id})
           x_init = [xy; v; b];
           P_init = eye(4,4);  
           agents{veh_id} = CoopAgent(veh_id, x_init, P_init, t, n_veh, system_setting);
           
            init_i = 0; % apply intialzation
            init_rewind = t;
        end        
            
        % Build prediction: F, f, Q
        agents{veh_id}.build_predict();         
       
        % Measurement\observation        
        if is_init == 1
%            agents{veh_id}.build_int_update('GPS', upt); % init mode
            agents{veh_id}.build_int_update('ZUPT', upt);
        else       
            is_gps = ~isempty(find(gps_agent == veh_id));
            if and(mod(iter_i, 4) == 0, is_gps)
                agents{veh_id}.build_int_update('GPS+IMU', upt);
            else
                agents{veh_id}.build_int_update('IMU', upt);            
            end
        end
                
    end    
    
    %% External observations: relative ranging       
     
    for i = 1 : length(veh_ids)
        
        veh_id = veh_ids(i);        
        xy_gt = agents{veh_id}.gt(end, :);
        
        fprintf('%i -> ', veh_id);
        
        % Add relative ranging
        if simulation_scenario == 3
            dist = sqrt( (epoch(:, 3) - xy_gt(1)).^2 + (epoch(:, 4) - xy_gt(2)).^2 ) ; 
            [vals, idx] = sort(dist);
            links_loc = epoch(idx(2:4), :);                                    % neighbors            
        else
            nidx = find( and( sqrt( (epoch(:, 3) - xy_gt(1)).^2 + (epoch(:, 4) - xy_gt(2)).^2 ) < system_setting.com_radius, epoch(:, 2) ~= veh_id )  ); 
            links_loc = epoch(nidx, :);                                         % neighbors        
        end

        for j = 1 : size(links_loc, 1)
            
            % Get indeces of neighbor
            veh_id_j = links_loc(j, 2);
            
            % This means that the vehicle is out of AOI
            if isempty(agents{veh_id_j})
                continue;
            end                      
            
            %agents{veh_id}.comm(agents{veh_id_j}, 'share-states');
            agents{veh_id}.comm(agents{veh_id_j}, 'UWB');            
            fprintf('%i ', veh_id_j);                        
        end        
        fprintf('\n');
        
        % Add infrastructure ranging
        if simulation_scenario ~= 3
            if ~isempty(infra_nodes)
                nidx = find( sqrt( (infra_nodes(:, 2) - xy_gt(1)).^2 + (infra_nodes(:, 3) - xy_gt(2)).^2 ) < system_setting.com_radius ); 
                for j = 1 : length(nidx)
                    agents{veh_id}.comm(infra_nodes(nidx(j), :), 'v2i');
                end
            end
        end
        
    end       
    fprintf('\n');
  
     if simulation_scenario == 3
         for i = 1 : size(infra_nodes, 1)
             dist = sqrt( (infra_nodes(i, 2) - epoch(:, 3)).^2 + (infra_nodes(i, 3) - epoch(:, 4)).^2 );
             [~, nidx] = min(dist);
             agents{epoch(nidx, 2)}.comm(infra_nodes(i, :), 'v2i');             
         end
     end
  
    %% Local filters
    for i = 1 : length(veh_ids)    

        veh_id = veh_ids(i);
        agent = agents{veh_id};
        if isempty(agent), continue; end
        
        % Assemble external and internal obs.
        H = [agent.H_int; agent.H_ext]; 
        z = [agent.z_int; agent.z_ext]; 
        R_add = diag( ones(length(agent.z_ext), 1)*system_setting.sigma_UWB^2 );
        R_int = agent.R_int;
        R = [R_int, zeros(size(R_int, 1), size(R_add, 2));  zeros(size(R_add, 1), size(R_int, 1)), R_add];
        h = afun_concat(agent.h_int, agent.h_ext);
        
        % Extended Kalman filter
        [x, P] = ekf_predict(agent.f, agent.F, agent.x, agent.P, agent.Q);
        [x, P] = ekf_update(x, z, h, H, P, R);
    
        %agent.x = x;
        %agent.P = P;
        agent.apply_update(x, P);
        
    end
    
 
    %% Covariance intersection  
    %if mod(t, 0.1) == 0
    new_agents{250} = [];
    for i = 1 : length(agents)
        if ~isempty(agents{i})
            new_agents{i} = copy(agents{i});
        end
    end
    
    if 1
        for i = 1 : length(veh_ids)

            agent = agents{veh_ids(i)};        
            idxs  = agent.idxs(1:4);
            %idxs  = 1:length(agent.x);                      
            %links = [ones(length(veh_ids), 1)*agent.id veh_ids];            % Communicating with all nodes
            
            %links = get_local_conn(agents, veh_ids, agent.id);
            %links = [ones(length(links), 1)*agent.id links];
            
            links = agent.links;                                           % Communicating only with neighbors
            links = [veh_ids(i), veh_ids(i); links];                      % Put the vehicles itsself into the links
            links = links(links(:,2)<1000, :);
            
            if size(links, 1) < 2
                continue
            end
            
            Ps = {};
            Is = {};
            for j = 1 : size(links, 1)     
                agent2 = agents{links(j, 2)};
                %H =  eye(length(idxs)) - agent2.K*[agent2.H_int; agent2.H_ext];
                %H =  [agent2.H_int; agent2.H_ext];
                %H = (H'*H)
                Ps{j} = agent2.P(idxs, idxs);
                Is{j} = inv(Ps{j});    
            end

            %ci_obj_func(ones(length(veh_ids), 1), Ps)
            n       = length(Is);
            w0      = ones(n, 1) / n;
            Aeq     = ones(1, n);
            beq     = 1;
            lb      = zeros(n, 1);
            ub      = ones(n, 1);    
            opts    = optimoptions('fmincon','Display','off');
            w       = fmincon(@(x) trace(inv(ci_obj_func(x, Is))), w0, [], [], Aeq, beq, lb, ub, [], opts);   

            %% Collect x and P, and fuse them
            % Nice examples of CI: https://arxiv.org/pdf/1610.01045.pdf        
            n = length(idxs);
            I_f = zeros(n, n);
            x_f = zeros(n, 1);
            I_avg = zeros(n, n);
            for j = 1 : size(links, 1)               
                I_i = Is{j};        
                I_f = I_f + w(j) * I_i;
                x_f = x_f + w(j) * I_i * agents{links(j, 2)}.x(idxs);
                I_avg = I_avg + I_i;        
            end
            P_f = inv(I_f);
            %P_f = inv(ci_obj_func(w, Ps));
            x_f = P_f * x_f;

            I_avg = I_avg / length(w);
            P_avg = inv(I_avg);

            % Covariance visualization
            if agent.id == sel_veh_id            
                figure(3); clf; hold on; 
                midxs = agents{sel_veh_id}.idxs;
                %midxs = midxs(1:2);
                midxs = 1:2;
                for j = 1 : size(links, 1)               
                    if min(eig(Ps{j})) < 0
                        error(sprintf('P is not positive definite! Link: %i -> %i', links(j, 1), links(j, 2)));
                    else
                        error_ellipse(Ps{j}(midxs,midxs), [0 0], 'style', 'b');                        
                    end      
                end
                
%                 for j = 1 : length(veh_ids)               
%                     Psl = agents{veh_ids(j)}.P(midxs,midxs);
%                     if min(eig(Psl)) < 0
%                         error(sprintf('P is not positive definite! Link: %i -> %i', links(j, 1), links(j, 2)));
%                     else
%                         error_ellipse(Psl, [0 0], 'style', 'b');                        
%                     end      
%                 end
                
                error_ellipse(P_f(midxs,midxs), [0 0], 'style', 'r--');
                error_ellipse(P_avg(midxs,midxs), [0 0], 'style', 'g');
                axis equal;
                
                %H = [agent.H_int; agent.H_ext];
                %return
            end

            % Back propagate to save time            
            new_agents{veh_ids(i)} .x(idxs) = x_f;
            new_agents{veh_ids(i)} .P(idxs, idxs) = P_f;

            %agents{veh_ids(i)}.apply_update(x_f, P_f);
%             for j = 1 : size(links, 1)    
%                 lidxs = agent.idxs;
%                 new_agents{links(j, 2)}.x = x_f;
%                 new_agents{links(j, 2)}.P = P_f;            
%             end
            
%             for j = 1 : size(veh_ids, 1)        
%                 agents{veh_ids(j)}.x(idxs) = x_f;
%                 agents{veh_ids(j)}.P(idxs, idxs) = P_f;            
%             end

        end
    end
    agents = new_agents;   
   
    
    %% Visualization and analytics
    
    % Vehicles positions and communication links
    figure(1); clf; hold on;
    if is_init
        title(sprintf("init t= %.3f\n", t));
    else
        title(sprintf("t= %.3f\n", t));
    end
    plot(roads(:, 1), roads(:, 2), 'k.');  
    if ~isempty(infra_nodes)
        plot(infra_nodes(:, 2), infra_nodes(:, 3), 'r.', 'MarkerSize', 15);    
    end
    for i = 1 : length(agents)        
        if isempty(agents{i}), continue; end
        
        x_hist    = agents{i}.my_x_hist();
        x_loc     = x_hist(end, :);
        xy_loc_gt = agents{i}.gt(end, :);
        z_loc     = agents{i}.z_int_hist(end, :);
        z_loc     = z_loc(agents{i}.sidx:(agents{i}.sidx+1));

        plot(x_loc(1), x_loc(2), 'b.', 'MarkerSize', 15);
        plot(xy_loc_gt(1), xy_loc_gt(2), 'g.', 'MarkerSize', 15);  
        plot(z_loc(1), z_loc(2), 'r.', 'MarkerSize', 15); 

        if agents{i}.id == sel_veh_id
            text(xy_loc_gt(1), xy_loc_gt(2)+20, sprintf('%i',  agents{i}.id), 'BackgroundColor', 'w', 'Color', 'r');
        else
            text(xy_loc_gt(1), xy_loc_gt(2)+20, sprintf('%i',  agents{i}.id), 'BackgroundColor', 'w');
        end
        
        if ~isempty(agents{i}.links)
            links_agents = agents{i}.links(agents{i}.links(:,2) < 1000, :);
            for j = 1 : size(links_agents, 1)
                nidx = links_agents(j, 2);
                xy_loc_gt_j = agents{nidx}.gt(end, :);
                plot([xy_loc_gt(1), xy_loc_gt_j(1)], [xy_loc_gt(2), xy_loc_gt_j(2)], 'g-');
            end

            links_infra = agents{i}.links(agents{i}.links(:,2) >= 1000, :);
            for j = 1 : size(links_infra, 1)
                nidx = find( links_infra(j, 2) == infra_nodes(:,1) );
                xy_loc_gt_j = infra_nodes(nidx, 2:3);
                plot([xy_loc_gt(1), xy_loc_gt_j(1)], [xy_loc_gt(2), xy_loc_gt_j(2)], 'r.-');
            end
       end
    end

    grid on; axis equal;    
    xlim(system_setting.AOI(1,1:2)); ylim(system_setting.AOI(2,1:2));       
    set(gca, 'FontSize', 12);    
    xlabel('[m]'); ylabel('[m]');
    
    % Selected vehicle's internal states    
    agent = agents{sel_veh_id};
    if ~isempty(agent)
        
        x_hist = agent.my_x_hist();
        
        idxs = agent.idxs;        
        n = size(x_hist, 1);
        dxy = (x_hist(end,1:4)-agent.gt(end, 1:4))';

        figure(2); clf; hold on;                

        subplot(4, 1, 1); hold on; 
        plot(1:n, x_hist(:,1), 'b.-' );
        plot(1:n, agent.gt(:, 1), 'g.-' );
        title(sprintf('[%i] X= %.3f', sel_veh_id, dxy(1)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 2); hold on; 
        plot(1:n,x_hist(:,2), 'b.-' );
        plot(1:n, agent.gt(:, 2), 'g.-' );
        title(sprintf('Y= %.3f', dxy(2)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 3); hold on; 
        plot(1:n, x_hist(:,3), 'b.-' );
        plot(1:n, agent.gt(:, 3), 'g.-' );
        title(sprintf('v= %.3f', dxy(3)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 4); hold on; 
        plot(1:n, x_hist(:,4)/pi*180, 'b.-' );
        plot(1:n, agent.gt(:, 4)/pi*180, 'g.-' );
        title(sprintf('bearing= %.3f', dxy(4)/pi*180));  
        set(gca, 'FontSize', 10);       
    end
        
    
    % Show sparisty of P
%     figure(3);
%     show_sparisty_pattern(P, 3);
    
    % Compate trajectries to ground truth   
    fprintf('Time: %.1f\n', t);    
    dxy_all = [];
    for i = 1 : length(veh_ids)
        if ~isempty(agents{veh_ids(i)})
            x_hist =  agents{veh_ids(i)}.my_x_hist();
            %dxy = agents{veh_ids(i)}.gt(:, 1:2) - x_hist(:, 1:2);
            dxy = agents{veh_ids(i)}.gt(end, 1:2) - x_hist(end, 1:2);
            dxy = sqrt(dxy(:, 1).^2 + dxy(:, 2).^2);
            xy_rmse = sqrt(sum(dxy.^2) / length(dxy));
            xy_std = std(dxy);
            fprintf('#%i RMSE: [%.3f]  STD: [%.3f]\n', veh_ids(i), xy_rmse, xy_std);
            %dxy_all  = [dxy_all; dxy(end)];
            dxy_all  = [dxy_all; dxy(end)];
        end
    end
    fprintf('TOTAL MEAN: [%.3f] MEDIAN: [%.3f]\n', mean(dxy_all), median(dxy_all));
    %input('Press key: ', 's');
    pause(viz_frame_rate);
end

sol.agents = agents;
sol.system_setting = system_setting;
sol.problem = problem;
%save(sprintf('solution_ci_%i_%i', simulation_scenario, system_setting.sigma_GPS*10), 'sol');
save(sprintf('solution_ci2_%i', simulation_scenario), 'sol')


