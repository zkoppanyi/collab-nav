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
agents = new_agents;   