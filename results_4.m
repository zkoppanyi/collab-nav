clear variables

simulation_settings;

sims{1} = 'results\solution_centralized_ukf';
sims{2} = 'results\solution_tracking';
sims{3} = 'results\solution_state_sharing';
sims{4} = 'results\solution_global_ci';
sims{5} = 'results\solution_local_ci';
sims{6} = 'results\solution_consensus';

form_t = 100; to_t = 110;
simulation_scenario = 4;

figure(1); clf; hold on;
cols = 'kmgbcr';

for sim_i = 1 : length(sims)
    sim_name = sims{sim_i};


        fname = sprintf('%s_%i', sim_name, simulation_scenario);
        load(fname);
       
        err = [];
        for i = 1 : length(sol.agents)
            agent = sol.agents{i};
            if isempty(agent), continue; end

            dt = agent.system_setting.dt;
            t = (agent.t0:dt:(agent.t0+size(agent.gt, 1)*dt-dt))';
            if sim_i < 2
                idxs = [1 2];
            else
                idxs = agent.idxs(1:2);                
            end
            
            idx = (agent.gt(:, end) == 0);
            %idx = idx(1:size(agent.x_hist, 1));
            dxy = agent.gt(idx, 1:2) - agent.x_hist(idx, idxs); 
            t = agent.gt(idx, 5);                                 
            
            dxy = sqrt(dxy(:, 1).^2 + dxy(:, 2).^2);   
            err = [err; t, dxy];                       
        end

        % Sum up 
        stat = [];
        for t = form_t : dt : to_t
           err_t = err( abs(err(:, 1) - t) < dt/2, 2);    
%            if size(err_t, 1) > 10
%                continue;
%            end
           stat = [stat; t, mean(err_t), median(err_t), std(err_t), length(err_t)];
        end
        stat = stat(~isnan(stat(:, 2)), :);

        B = 1/9*ones(9,1);    
        err_filt = filter(B,1,stat(:,2));
        err_filt_std = filter(B,1,stat(:,4));
        
        plot(stat(:,1)-form_t, err_filt(:), [cols(sim_i), '-'], 'LineWidth', 1.5);        
        %plot(stat(idx,1)-form_t, err_filt(idx) + err_filt_std(idx) , [cols(sim_i), '--'], 'LineWidth', 0.5);
        %plot(stat(idx,1)-form_t, err_filt(idx) - err_filt_std(idx), [cols(sim_i), '--'], 'LineWidth', 0.5);

end

set(gca, 'FontSize', 14);
xlabel('Time [s]');
ylabel('Average RMSE [m]');
grid on
ylim([0 1])
ylabel('\Delta [m]');
plot([0 10], [0.3 0.3], 'k:', 'LineWidth', 1);
legend('Centralized', 'Tracking', 'State sharing', 'Global CI', 'Local CI', 'Consensus')




