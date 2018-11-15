
function [v, G, A] = get_local_conn(agents, veh_ids, v_start)

    edges = [];
    for i = 1:length(veh_ids)
        if ~isempty(agents{veh_ids(i)})
            links = agents{veh_ids(i)}.links;
            if ~isempty(links)
                links = links( links(:, 2) < 1000, :);
                for j = 1 : size(links, 1)
                    idx2 = find(links(j, 2) == veh_ids);
                    edges = [edges; i idx2];
                end
            end
        end
    end
    G = graph(edges(:,1), edges(:,2));
    v = bfsearch(G, v_start);
    A = adjacency(G);