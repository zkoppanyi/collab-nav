function L = get_laplacian(curr_id, agents, veh_ids, G)

    n_agents = length(veh_ids);
    A = zeros(n_agents, n_agents);
    for i = 1 : n_agents
        
        agent = agents{veh_ids(i)};
        if isempty(shortestpath(G, curr_id, agent.id ))
            continue;
        end
            
        links = agent.links;  
        if isempty(links), continue; end
        for j = 1 : size(links, 1)
            if links(j, 2) >= 1000,  continue; end
            idx = find(links(j, 2) == veh_ids);
            A(i, idx) = 1; A(idx, i) = 1;
        end
    end

    C = [];
    for i = 1 : n_agents
        for j = (i+1) : n_agents
            if A(i,j) == 1
                col = zeros(n_agents, 1);
                col(i) = 1;
                col(j) = -1;
                C = [C, col];
            end
        end
    end
    L = C*C';

    %W = A./sum(A, 2);
    %W*ones(size(W,1), 1);        