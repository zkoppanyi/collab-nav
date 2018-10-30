clear variables

simulation_settings;


sims{1} = 'solution_centralized';
sims{2} = 'solution_ci';
sims{3} = 'solution_consensus';


form_t = 100; to_t = 130;
simulation_scenario = 3;

figure(1); clf; hold on;
cols = 'kbgmr';

for sim_i = 1 : length(sims)
    sim_name = sims{sim_i};


        fname = sprintf('%s_%i', sim_name, simulation_scenario);
        load(fname);

        if (sim_i == 1)
            sol = solution_centralized;  
        end
        
        err = [];
        for i = 1 : length(sol.agents)
            agent = sol.agents{i};
            if isempty(agent), continue; end

            dt = agent.system_setting.dt;
            t = (agent.t0:dt:(agent.t0+size(agent.gt, 1)*dt-dt))';
            if sim_i == 1
                idxs = [1 2];
            else
                idxs = agent.idxs(1:2);                
            end
            
            if sim_i >= 2
                idx = (agent.gt(:, end) == 0);
                idx = idx(1:size(agent.x_hist, 1));
                dxy = agent.gt(idx, 1:2) - agent.x_hist(idx, idxs); 
                t = agent.gt(idx, 5);                                 
            else            
                dxy = agent.gt(:, 1:2) - agent.x_hist(:, idxs);                                  
            end
            
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

        idx = find(stat(:,2) < Inf);
        %plot(stat(idx,1), stat(idx,3), [cols(test_i) '.-']);
        n_veh = medfilt1(stat(:,4), 5);

        B = 1/9*ones(9,1);    
        err_filt = filter(B,1,stat(:,3));
        %err_filt = medfilt1(stat(:,3), 5);
        %err_filt = stat(:,2);
        
        plot(stat(idx,1)-form_t, err_filt(idx), [cols(sim_i), '-'], 'LineWidth', 1);
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
ylim([0 2])
ylabel('\Delta [m]');
legend('Centralized', 'State sharing', 'Centralized CI', 'Scalable CI', 'Consensus')




