clear variables

simulation_settings;

sims{1} = 'results\solution_centralized';
sims{2} = 'results\solution_single';
sims{3} = 'results\solution_tracking';
sims{4} = 'results\solution_state_sharing';
sims{5} = 'results\solution_global_ci';
sims{6} = 'results\solution_local_ci';
sims{7} = 'results\solution_consensus';

form_t = 100; to_t = 110;
simulation_scenario = 2;

figure(1); clf; hold on;
cols = 'kymgbcr';

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
            if sim_i < 3
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
           stat = [stat; t, mean(err_t), median(err_t), length(err_t)];
        end
        stat = stat(~isnan(stat(:, 2)), :);

        idx = find(stat(:,2) < Inf);
        %plot(stat(idx,1), stat(idx,3), [cols(test_i) '.-']);
        n_veh = medfilt1(stat(:,4), 5);

        B = 1/9*ones(9,1);    
        err_filt = filter(B,1,stat(:,2));
        %err_filt = medfilt1(stat(:,3), 5);
        %err_filt = stat(:,2);
        
        plot(stat(idx,1)-form_t, err_filt(idx), [cols(sim_i), '-'], 'LineWidth', 1.5);
        %[hAx,hLine1,hLine2] = plotyy(stat(idx,1), err_filt, stat(:,1), n_veh);
        %hLine1.LineStyle = '-';
        %hLine1.Color = cols(test_i);
        %hLine2.LineStyle = '-';
        %hLine2.Color = 'k';

end

set(gca, 'FontSize', 14);
xlabel('Time [s]');
ylabel('Average RMSE [m]');
grid on
ylim([0 5])
ylabel('\Delta [m]');
plot([0 10], [5.0 5.0], 'k--', 'LineWidth', 0.5);
plot([0 10], [0.3 0.3], 'k:', 'LineWidth', 1);
legend('Centralized', 'No coop.', 'Tracking', 'State sharing', 'Global CI', 'Local CI', 'Consensus')




