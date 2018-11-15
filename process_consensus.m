% Consensus filter
% Code based on:
% https://ac.els-cdn.com/S0005109816300188/1-s2.0-S0005109816300188-main.pdf?_tid=d446a434-ef8e-48ca-a025-be3f6977a7a4&acdnat=1539389828_3e199ac4f616bfd57f1dbf569fdb3b61    

[~, G, Adj] = get_local_conn(agents, veh_ids, veh_ids(1)); 
%C = create_consensus_matrix(Adj, 'laplacian-const');
C = create_consensus_matrix(Adj, 'max-degree');
%C = create_consensus_matrix(Adj, 'laplacian-vary');

for i = 1 : length(veh_ids)    

    agent = agents{veh_ids(i)};
    agent.build_predict(); 

    [h, H, z, R] = agent.build_update();                    
    z = z - afun_eval(h, agent.x) + H*agent.x;
    agent.dq     = H'*inv(R)*z;
    agent.dOmega = H'*inv(R)*H;

end 

% Simple average to test consensus
%     w = 1/length(veh_ids);
%     q_sum       = w * agents{veh_ids(1)}.q;
%     Omega_sum   = w * agents{veh_ids(1)}.Omega;
%     dq_sum      = w * agents{veh_ids(1)}.dq;
%     dOmega_sum  = w * agents{veh_ids(1)}.dOmega;
%     for i = 2 : length(veh_ids)
%         q_sum       = q_sum      + w * agents{veh_ids(i)}.q;
%         Omega_sum   = Omega_sum  + w * agents{veh_ids(i)}.Omega;
%         dq_sum      = dq_sum     + w * agents{veh_ids(i)}.dq;
%         dOmega_sum  = dOmega_sum + w * agents{veh_ids(i)}.dOmega;
%     end
%     for i = 1 : length(veh_ids)
%         agents{veh_ids(i)}.q        = q_sum;
%         agents{veh_ids(i)}.Omega    = Omega_sum;
%         agents{veh_ids(i)}.dq       = dq_sum;
%         agents{veh_ids(i)}.dOmega   = dOmega_sum;
%     end

% Consensus
Omegas_est = {};

iter = [];
%for ci = 1 : length(C)
for ci = 1 : 50

     new_agents = agents;
     iter_col = zeros(length(veh_ids), 1);

     for i = 1 : length(veh_ids)

        agent = agents{veh_ids(i)};
        links = agent.links;          
        if isempty(links), continue; end            
        
        links = links(links(:, 2) < 1000, :);
        if isempty(links), continue; end            

        % Save iterations
        idxs = agents{sel_veh_id}.idxs;
        iter_col(i)     = agent.q(idxs(1));
        iter_col(i)     = norm(agent.q);

        %w = C{ci}(i,i);
        w = C(i,i); 
        sum_w = w;
        q      = w * agent.q;
        Omega  = w * agent.Omega;            
        dq     = w * agent.dq;
        dOmega = w * agent.dOmega;            

        for j = 1 : size(links, 1)
            idx = find(links(j, 2) == veh_ids);
            agent2 = agents{links(j, 2)};

            %w = C{ci}(i, idx);
            w      = C(i, idx); 
            sum_w = sum_w + w;
            q      = q       + w * agent2.q;
            Omega  = Omega   + w * agent2.Omega;
            dq     = dq      + w * agent2.dq;
            dOmega = dOmega  + w * agent2.dOmega;    
        end

        agent.q         = q;
        agent.Omega     = Omega;
        agent.dq        = dq;
        agent.dOmega    = dOmega;         

        if sum(isnan(q)) ~= 0
            error('Somthing went wrong! NaN value in q!')
        end

        if abs(sum_w - 1) > 1e-5
            error('Somthing went wrong! sum_w ~= 1!')                
        end

        Omegas_est{ci} = inv(Omega);            
        agents{veh_ids(i)} = agent;        

     end         

     iter = [iter, iter_col];
     agents = new_agents;                     
end

figure(3); clf; hold on;
plot(iter');
title ('Convergence');


%gamma = length(veh_ids);
gamma = 1;
Omegas = {};
for i = 1 : length(veh_ids)    

    agent = agents{veh_ids(i)};

    agent.q     = agent.q       + gamma*agent.dq;
    agent.Omega = agent.Omega   + gamma*agent.dOmega;
    P = inv(agent.Omega);
    x = P*agent.q;

    agent.apply_update(x, P);
    agent.build_predict();  

    x_hat = afun_eval(agent.f, x);
    A = agent.F;
    W = inv(agent.Q);
    agent.Omega = W - W*A*inv(agent.Omega + A'*W*A)*A'*W;
    agent.q = agent.Omega*x_hat;
    Omegas{i} = inv(agent.Omega);        


    agents{veh_ids(i)} = agent;        
end 

%% Covariances
%     figure(3); clf; hold on; 
%     sidx = agents{sel_veh_id}.sidx;
%     for i = 1 : length(Omegas)
%         error_ellipse(Omegas{i}(sidx:(sidx+1), sidx:(sidx+1)), [0 0], 'style', 'b');
%         error_ellipse(Omegas_est{i}(sidx:(sidx+1), sidx:(sidx+1)), [0 0], 'style', 'r');
%     end
%     axis equal;