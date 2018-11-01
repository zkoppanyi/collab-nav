clear variables

simulation_settings;

%gps_error_tests = [0.1 0.5 1.0 5.0 10.0];
gps_error_tests = 5.0;

sims{1} = 'results\solution_centralized';
sims{2} = 'results\solution_single';
sims{3} = 'results\solution_tracking';
sims{4} = 'results\solution_state_sharing';
sims{5} = 'results\solution_ci_stat';
sims{6} = 'results\solution_consensus_stat';

form_t = 100; to_t = 110;
simulation_scenario = 1;

figure(1); clf; hold on;
plot([0 to_t-form_t], [gps_error_tests gps_error_tests], 'k--', 'LineWidth', 1.5)

cols = 'kmybgr';
lw{1} = '-';
lw{2} = '-';
lw{3} = '-';
lw{4} = '-';
lw{5} = '-';
lw{6} = '-';
for sim_i = 1 : length(sims)
    sim_name = sims{sim_i};
    for test_i = 1 : length(gps_error_tests)
        fname = sprintf('%s_%i_%i', sim_name, simulation_scenario, gps_error_tests(test_i)*10);
        load(fname);

        if (sim_i == 1)
            sol = solution_centralized;            
        elseif (sim_i == 2)
            sol = solution_single;                 
        elseif (sim_i == 3)
            sol = solution_state_sharing;  
        elseif (sim_i == 4)
            sol = solution_state_sharing;          
        end
        
        err = [];
        for i = 1 : length(sol.agents)
            agent = sol.agents{i};
            if isempty(agent), continue; end

            dt = agent.system_setting.dt;
            t = (agent.t0:dt:(agent.t0+size(agent.gt, 1)*dt-dt))';
            if or((sim_i == 1), (sim_i == 2))
                idxs = [1 2];
            else
                idxs = agent.idxs(1:2);
            end
            
            if or( sim_i == 5, sim_i == 6)
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
           if isnan(mean(err_t))
                continue
           end
           stat = [stat; t, mean(err_t), median(err_t), length(err_t)];
        end

        idx = find(stat(:,2) < Inf);
        %plot(stat(idx,1), stat(idx,3), [cols(test_i) '.-']);
        n_veh = medfilt1(stat(:,4), 5);

        B = 1/10*ones(10,1);    
        err_filt = filter(B,1,stat(:,2));
        %err_filt = medfilt1(stat(:,3), 5);
        %err_filt = stat(:,2);

        plot(stat(idx,1)-form_t, err_filt(idx), [cols(sim_i), lw{test_i}], 'LineWidth', 1.5);
        %[hAx,hLine1,hLine2] = plotyy(stat(idx,1), err_filt, stat(:,1), n_veh);
        %hLine1.LineStyle = '-';
        %hLine1.Color = cols(test_i);
        %hLine2.LineStyle = '-';
        %hLine2.Color = 'k';
        
    end
end

set(gca, 'FontSize', 14);
xlabel('Time [s]');
ylabel('\Delta [m]');
grid on
legend('GPS error', 'Centralized', 'No coop.', 'Tracking', 'State sharing', 'CI', 'Consensus EKF');
ylim([0 6])



