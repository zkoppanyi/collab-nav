
function [v, G] = get_local_conn(agents, veh_ids, v_start)

    edges = [];
    for i = 1:length(veh_ids)
        if ~isempty(agents{veh_ids(i)})
            edges = [edges; agents{veh_ids(i)}.id, agents{veh_ids(i)}.id; agents{veh_ids(i)}.links];
        end
    end
    G = graph(edges(:,1), edges(:,2));
    v = bfsearch(G, v_start);