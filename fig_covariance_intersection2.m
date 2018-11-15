clear variables

simulation_settings;

sims{1} = 'results\solution_centralized';
%sims{2} = 'results\solution_single';
sims{2} = 'results\solution_tracking';
sims{3} = 'results\solution_state_sharing';
sims{4} = 'results\solution_global_ci';
%sims{5} = 'results\solution_local_ci';
sims{5} = 'results\solution_consensus';

t = 104.0;
simulation_scenario = 4;

figure(1); clf; hold on;
cols = 'kmgbr';

for sim_i = 1 : length(sims)
    sim_name = sims{sim_i};


        fname = sprintf('%s_%i', sim_name, simulation_scenario);
        load(fname);

        agent = sol.agents{6};
        idx = find( abs(agent.gt(:, 5) - t) < 0.05/2 );
        
        if isempty(idx)
            disp('No data!');
            continue;
        end
        
        if sim_i == 1
            idxs = [1 2];
            P = agent.x_hist(idx, 5:end);
            P = reshape(P, 4, 4);
        else
            idxs = agent.idxs(1:2);
            P = agent.x_hist(idx, (7*4+1):end);            
            n = sqrt(length(P));
            P = reshape(P, n, n);
        end
        P = P(idxs, idxs);
        
        if min(eig(P)) > 0        
            error_ellipse(P, [0 0], 'style', cols(sim_i)); 
        end
        
end
legend('Centralized', 'Tracking', 'State sharing', 'Global CI', 'Consensus')

axis equal;
grid on;
set(gca, 'FontSize', 14)






