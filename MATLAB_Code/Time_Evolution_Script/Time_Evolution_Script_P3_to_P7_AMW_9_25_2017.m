% This code tests under what conditions, if any, it is possible to reach a
% condition similar to that observed at P7 beginning at the connectivity
% observed at P3.

% For each time step, will 1. remove synapses stochastically.

% After this step, axons will be checked and any axon whose number of
% synapses with the seed cell (PC1, which is represented by column 1 in
% the connectivity matrices we will be working with) goes to zero will be
% removed from the connectivity matrix, so the results reflect what we
% would have measured by finding inputs connected only to the seed cell at
% the simulated "P7" timepoint.

% 2. Add synapses stochastically

% Added synapses will be added ONLY to connections that already exist
% between axons and PCs.

% An implication of this order of operations is that once an axon becomes
% disconnected from a target cel, it cannot re-connect with that target
% cell under these conditions.

% At each step the connectivity matrix will be compared to the observed P7
% connectivity matrix by several criteria.

% The simulation will end when either the criteria converge with one
% another to within specified errors OR the maximum number of timesteps has
% been reached.

% NOTE: "power_targ_curr" is the same as gamma.

% This script makes use of "createPatches.m"
% (https://www.mathworks.com/matlabcentral/answers/uploaded_files/30414/createPatches.m)

% Alyssa Wilson, October 2016
% Lichtman Lab, Harvard University

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load observed P3 and P7 connectivity matrices

% Load P3 data
% File is "P3_Observed_PC_Connectivity_Synapse_Numbers.mat"
% Variable is "P3_PCconnectivity". Rows are axons, columns are PCs.
load 'I:\P3_P7_Comparisons\Time_Evolution_Script_for_Pub\P3_Observed_PC_Connectivity_Synapse_Numbers.mat';

% Load P7 data
% File is "P7_Observed_PC_Connectivity_Synapse_Numbers.mat"
% Variable is "P7_PCconnectivity". Rows are axons, columns are PCs.
load 'I:\P3_P7_Comparisons\Time_Evolution_Script_for_Pub\P7_Observed_PC_Connectivity_Synapse_Numbers.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

N_ITER_PER_PARAM_SET         = 10;

P_REMOVAL_TO_TEST            = 0.01; %[0 0.0005 0.001 0.005 0.01];
POWER_TARGET_PROB            = 1.3; %[0 1.0:0.1:2.0 2.5 3.0 3.5 4.0];

alltrialshasconverged        = [];
alltrialsnstepstoconv        = [];

alltrialspremoval            = [];
alltrialspowtarprob          = [];

groupedrepeatiter_premoval   = [];
groupedrepeatiter_powtarprob = [];
groupedrepeatiter_meanhascon = [];
groupedrepeatiter_meantcon   = [];

% Columns for this array are iteration number, cf number, average slope of
% number of synapses total as a function of time step
avgrateofsynchangepercfperiter = [];


count = 1;
for tparam1 = 1:length(P_REMOVAL_TO_TEST)
    p_rem_curr  = P_REMOVAL_TO_TEST(tparam1);
    
    for tparam2 = 1:length(POWER_TARGET_PROB) % FOR SPECIAL TESTS OF PARTICULAR PARAMETER PAIRS, CHANGE UPPER LIMIT TO 1
        power_targ_curr = POWER_TARGET_PROB(tparam2); % FOR SPECIAL TESTS OF PARTICULAR PARAMETER PAIRS, CHANGE INDEX TO tparam1
        
        hasconvforparam = [];
        tconvforparam   = [];
        trackerrs       = [];
        
        disp(sprintf('Running new parameter set %d',count));
        count           = count + 1;
        
        % Set up array to hold quantiles for the time-evolved
        % connectivity matrix
        % Each row is an iteration
        % Each column is a quantile: 0.25 0.50 0.75 0.95
        quantilesalliter = [];
        cumstocalculate  = [0.25 0.50 0.75 0.95];
        
        % Set up an array to hold the quantiles for the observed P7 data
        quantilesP7      = [];
        
        for itercurr = 1:N_ITER_PER_PARAM_SET
            disp(sprintf('Running iteration %d',itercurr));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % Set up arrays that track the addition of synapses over time
            
            % Set up an array to track the number of synapses added per
            % timestep
            nsyns_added_per_timestep = [];
            
            % An array to keep track of the total number of synapses made
            % by each climbing fiber over timesteps
            % Columns are the climbing fiber unique ID, the timestep
            % number, and the total number of synapses made by that
            % climbing fiber at the current timestep
            nsyns_tot_for_cfs_over_timesteps = [];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set up initial conditions by calculating comparison values from P3
            % connectivity matrix
            
            connectivity_curr          = P3_PCconnectivity;
            
            permrowids                 = [1:size(P3_PCconnectivity,1)]'; % Make a column vector
            permcolids                 = [1:size(P3_PCconnectivity,2)]'; % Make a column vector
            
            % Calculate median and skewness about mean of histogram
            allPCconnectivitynozeros_i = P3_PCconnectivity(:);
            allPCconnectivitynozeros_i(find(~allPCconnectivitynozeros_i)) = [];
            hAllConnectivity_i         = histogram(allPCconnectivitynozeros_i);
            
            ninputs_P3_i               = size(P3_PCconnectivity,1);
            
            nPCs_P3_i                  = size(P3_PCconnectivity,2);
            
            median_P3_i                = median(allPCconnectivitynozeros_i);
            
            mean_P3_i                  = mean(allPCconnectivitynozeros_i);
            stdev_P3_i                 = std(allPCconnectivitynozeros_i); % Used for standardizing the measurement of the skewness
            skewness_P3_i              = mean( ((allPCconnectivitynozeros_i - mean_P3_i)./stdev_P3_i).^3 );
            
            % Calculate mean and standard deviation of total numbers of synapses made
            % per axon
            nsynstot_P3_i              = sum(P3_PCconnectivity,2);
            mean_nsynstot_P3_i         = mean(nsynstot_P3_i);
            stddev_nsynstot_P3_i       = std(nsynstot_P3_i);

                        
            % Calculate number of cells contacted by each axon
            binary_P3 = P3_PCconnectivity ~= 0;
            n_cells_contacted_per_axon_i   = sum(binary_P3,2);
            mean_ncells_P3_i               = mean(n_cells_contacted_per_axon_i);
            std_ncells_P3_i                = std(n_cells_contacted_per_axon_i);

            
            % Calculate number of axons contacting each cell
            n_axons_per_cell_i             = sum(binary_P3,1);
            mean_nax_P3_i                  = mean(n_axons_per_cell_i);
            std_nax_P3_i                   = std(n_axons_per_cell_i);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate values for observed P7 connectivity matrix that will be used as
            % a standard to compare to
            
            % Calculate median and skewness about mean of histogram
            allPCconnectivitynozeros_P7 = P7_PCconnectivity(:);
            allPCconnectivitynozeros_P7(find(~allPCconnectivitynozeros_P7)) = [];
            hAllConnectivity_P7         = histogram(allPCconnectivitynozeros_P7);
            
            quantilesP7                 = quantile(allPCconnectivitynozeros_P7,cumstocalculate);
            
            ninputs_P7                  = size(P7_PCconnectivity,1);
            
            nPCs_P7                     = size(P7_PCconnectivity,2);
            
            median_P7                   = median(allPCconnectivitynozeros_P7);
            
            mean_P7                     = mean(allPCconnectivitynozeros_P7);
            stdev_P7                    = std(allPCconnectivitynozeros_P7);
            skewness_P7                 = mean( ((allPCconnectivitynozeros_P7 - mean_P7)./stdev_P7).^3 );
            
            % Calculate mean and standard deviation of total numbers of synapses made
            % per axon
            nsynstot_P7                 = sum(P7_PCconnectivity,2);
            mean_nsynstot_P7            = mean(nsynstot_P7);
            stddev_nsynstot_P7          = std(nsynstot_P7);

            
            % Calculate number of cells contacted by each axon
            binary_P7 = P7_PCconnectivity ~= 0;
            n_cells_contacted_per_axon_P7  = sum(binary_P7,2);
            mean_ncells_P7                 = mean(n_cells_contacted_per_axon_P7);
            std_ncells_P7                  = std(n_cells_contacted_per_axon_P7);
            
            % Calculate number of axons contacting each cell
            n_axons_per_cell_P7            = sum(binary_P7,1);
            mean_nax_P7                    = mean(n_axons_per_cell_P7);
            std_nax_P7                     = std(n_axons_per_cell_P7);
            
            close all;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            % Define minimum errors or amounts of overlap that should all be reached before simulation is
            % considered to have converged to an observed-P7-like condition
            epsilon_ninputs     = 20.0;
            epsilon_nPCs        = 12.0;
            
            epsilon_median      = 1.0;
            epsilon_skewness    = 1.0;
            
            epsilon_mean_ntot   = 4.0; 
            epsilon_std_ntot    = 3.0;

            epsilon_mean_ncells = 1.0;
            epsilon_std_ncells  = 1.0;
        
            epsilon_mean_nax    = 7;
            epsilon_std_nax     = 9;
            
            % Initialize current errors between initial conditions and observed P7
            % conditions
            err_ninputs_curr    = abs(ninputs_P7 - ninputs_P3_i);
            err_nPCs_curr       = abs(nPCs_P7 - nPCs_P3_i);
            err_median_curr     = abs(median_P7 - median_P3_i);
            err_skewness_curr   = abs(skewness_P7 - skewness_P3_i);

            err_mean_ntot_curr   = abs(mean_nsynstot_P7 - mean_nsynstot_P3_i);
            err_std_ntot_curr    = abs(stddev_nsynstot_P7 - stddev_nsynstot_P3_i);

            err_mean_ncells_curr = abs(mean_ncells_P7 - mean_ncells_P3_i);
            err_std_ncells_curr  = abs(std_ncells_P7 - std_ncells_P3_i);

            err_mean_nax_curr    = abs(mean_nax_P7 - mean_nax_P3_i);
            err_std_nax_curr     = abs(std_nax_P7 - std_nax_P3_i);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set up timestepping
            
            t       = 1; % index of timestep
            tmax    = 5000; % Maximum number of timesteps that will be taken before loop is stopped
            minconv = 1220; % Minimum timestep for convergence used when fine-tuning to find
            % the timestep where best accuracy is achieved (based on initial convergence results).
            maxconv = 1380; % Maximum timestep for timestep fine-tuning (based on initial convergence results).
            hasconverged = 0;
            
            prob_remove_syn = p_rem_curr;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Start time stepping
            
            while (t<=tmax && ~hasconverged) % FOR FINE-TUNING, USE (t<=tmax) && (~hasconverged || (t<minconv || t>maxconv) )
                % Update the tracking array with the total number of
                % synapses made per climbing fiber in the current timestep
                nsyns_tot_for_cfs_over_timesteps = [nsyns_tot_for_cfs_over_timesteps; [permrowids repmat(t,size(permrowids)) sum(connectivity_curr,2)]];
                
                % Update the array that tracks the errors for a simulation.
                % Each column is the timestepping behavior for a single error
                % parameter.
                trackerrs = [trackerrs; [err_ninputs_curr err_nPCs_curr err_median_curr err_skewness_curr err_mean_ntot_curr ...
                    err_std_ntot_curr err_mean_ncells_curr err_std_ncells_curr err_mean_nax_curr err_std_nax_curr]];
                
                
                % Get indices of current non-zero elements in the connectivity matrix
                connected_el_ids   = find(connectivity_curr);
                
                % Determine whether to remove any synapses from any connected elements
                selection_criteria_decrement = rand(size(connected_el_ids));
                
                weights_for_per_syn_removal  = repmat(sum(connectivity_curr,2),[1, size(connectivity_curr,2)]);
                weights_for_per_syn_removal  = weights_for_per_syn_removal./sum(sum(connectivity_curr));
                weights_selection_criteria   = weights_for_per_syn_removal(connected_el_ids);
                
                % Combine each random number with the correct weight for
                % the connection you're drawing it for
                selection_criteria_decrement = selection_criteria_decrement./weights_selection_criteria;
                
                % If the value of any of these elements is less than the probability of
                % a removal of a synapse, mark that as a connection to be decremented
                % by one synapse
                idsfordecrement    = connected_el_ids(find(selection_criteria_decrement <= prob_remove_syn));

                % Decrement the specified elements of the current connectivity matrix
                % by one synapse
                connectivity_curr(idsfordecrement) = connectivity_curr(idsfordecrement) - 1;
                
                
                % Check for zeros in the seed cell column and remove any rows with
                % zeros here
                disconnectedrows = find(connectivity_curr(:,1) == 0);
                connectivity_curr(disconnectedrows,:) = [];
                permrowids(disconnectedrows) = [];
                
                % Remove any columns in the connectivity matrix that are all zeros
                disconnectedcolumns = find(sum(connectivity_curr,1)==0);
                connectivity_curr(:,disconnectedcolumns) = [];
                permcolids(disconnectedcolumns) = [];

                % Get indices of remaining non-zero elements in connectivity matrix
                connected_el_ids = find(connectivity_curr);
                
                
                % Choose which axons get synapses added using a weighted probability
                % that is equivalent to the fraction of synapses in total that are made
                % by each axon
                nsynsperaxtotcurr                  = sum(connectivity_curr,2);
                probsforaxonscurr                  = nsynsperaxtotcurr./(sum(nsynsperaxtotcurr));

                selection_criteria_increment       = rand(size(connectivity_curr,1),1);
                axsforincrement                    = find(selection_criteria_increment <= probsforaxonscurr);
                
                % Calculate the total number of synapses added in the
                % current time step
                nsyns_added_per_timestep           = [nsyns_added_per_timestep; length(axsforincrement)];
                
                % For the axons that are chosen for an increment, decide which
                % connection gets the extra synapse. Probabilities will be weighted by
                % the relative number of synapses each connection has
                for j = 1:length(axsforincrement)
                    axrowcurr     = axsforincrement(j);
                    colids        = [1:length(connectivity_curr(axrowcurr,:))];
                    probscurrnum  = connectivity_curr(axrowcurr,:);
                    probscurrnum  = probscurrnum.^(power_targ_curr);
                    totsynscurr   = sum(probscurrnum);
                    probscurr     = probscurrnum./totsynscurr;
                    zeroscurr     = find(~probscurr);
                    probscurr(zeroscurr) = [];
                    colids(zeroscurr)    = [];
                    
                    selector        = rand(1);
                    idinterval      = 1;
                    cumprobcurr     = probscurr(idinterval);
                    
                    while selector > cumprobcurr
                        idinterval   = idinterval + 1;
                        cumprobcurr  = cumprobcurr + probscurr(idinterval);
                    end
                    
                    conntoaddtocurr = colids(idinterval);
                    connectivity_curr(axrowcurr,conntoaddtocurr) = connectivity_curr(axrowcurr,conntoaddtocurr) + 1;
                end
                
                
                % Calculate test quantities of new connectivity matrix
                ninputs_curr                        = size(connectivity_curr,1);
                nPCs_curr                           = size(connectivity_curr,2);
                
                connectivity_curr_nozeros           = connectivity_curr(:);
                connectivity_curr_nozeros(find(~connectivity_curr_nozeros)) = [];
                median_P3_curr                      = median(connectivity_curr_nozeros);
                stdev_P3_curr                       = std(connectivity_curr_nozeros);
                skewness_P3_curr                    = mean( ((connectivity_curr_nozeros - mean(connectivity_curr_nozeros))./stdev_P3_curr).^3 );
                
                nsynstot_P3_curr                    = sum(connectivity_curr,2);
                mean_nsynstot_P3_curr               = mean(nsynstot_P3_curr);
                stddev_nsynstot_P3_curr             = std(nsynstot_P3_curr);
                
                binary_P3_curr                      = connectivity_curr ~= 0;
                n_cells_contacted_per_axon_P3_curr  = sum(binary_P3_curr,2);
                mean_ncells_P3_curr                 = mean(n_cells_contacted_per_axon_P3_curr);
                std_ncells_P3_curr                  = std(n_cells_contacted_per_axon_P3_curr);
                
                n_axons_per_PC_P3_curr              = sum(binary_P3_curr,1);
                mean_nax_P3_curr                    = mean(n_axons_per_PC_P3_curr);
                std_nax_P3_curr                     = std(n_axons_per_PC_P3_curr);
                
                % Calculate errors
                err_ninputs_curr     = abs(ninputs_P7 - ninputs_curr);
                err_nPCs_curr        = abs(nPCs_P7 - nPCs_curr);
                err_median_curr      = abs(median_P7 - median_P3_curr);
                err_skewness_curr    = abs(skewness_P7 - skewness_P3_curr);

                err_mean_ntot_curr   = abs(mean_nsynstot_P7 - mean_nsynstot_P3_curr);
                err_std_ntot_curr    = abs(stddev_nsynstot_P7 - stddev_nsynstot_P3_curr);

                err_mean_ncells_curr = abs(mean_ncells_P7 - mean_ncells_P3_curr);
                err_std_ncells_curr  = abs(std_ncells_P7 - std_ncells_P3_curr);

                err_mean_nax_curr    = abs(mean_nax_P7 - mean_nax_P3_curr);
                err_std_nax_curr     = abs(std_nax_P7 - std_nax_P3_curr);

                
                % Check whether all errors are less than cutoffs, and if so, set
                % hasconverged to 1 and stop timestepping
                if (err_ninputs_curr<=epsilon_ninputs && err_nPCs_curr<=epsilon_nPCs ...
                        && err_median_curr<=epsilon_median ...
                        && err_skewness_curr<=epsilon_skewness ...
                        && err_mean_ntot_curr<=epsilon_mean_ntot && err_std_ntot_curr<=epsilon_std_ntot...
                        && err_mean_ncells_curr<=epsilon_mean_ncells && err_std_ncells_curr<=epsilon_std_ncells ...
                        && err_mean_nax_curr<=epsilon_mean_nax &&err_std_nax_curr<=epsilon_std_nax)

                    hasconverged = 1;
                    disp(sprintf('Simulation has converged at timestep %d.',t));
                    %t=t+1; % UNCOMMENT THIS LINE FOR FINE-TUNING
                elseif (isempty(connectivity_curr))
                    % If all the entries for the connectivity matrix have been removed,
                    % break the loop
                    t=tmax + 2; % Add 2 instead of one here so you can tell (in principle) when this is the case
                else
                    % Continue onto next timestep (only do this if solution has not
                    % converged, so you can tell at exactly which timestep the solution
                    % converges when it does)
                    % hasconverged = 0; % UNCOMMENT THIS LINE FOR
                    % FINE-TUNING
                    t=t+1;
                end
                
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            alltrialshasconverged     = [alltrialshasconverged; hasconverged];
            alltrialsnstepstoconv     = [alltrialsnstepstoconv; t];
            
            alltrialspremoval         = [alltrialspremoval; p_rem_curr];
            alltrialspowtarprob       = [alltrialspowtarprob; power_targ_curr];
            
            hasconvforparam           = [hasconvforparam; hasconverged];
            tconvforparam             = [tconvforparam; t];
                        
            
            % Compute the quantiles for only runs that have converged
            if hasconverged
                quantilesalliter          = [quantilesalliter; quantile(connectivity_curr_nozeros,cumstocalculate)];
            end
            
            % Plot the total number of synapses made by each climbing fiber over time
            cfstoplotinfofor   = unique(nsyns_tot_for_cfs_over_timesteps(:,1));
            figure()
            for plottingid = 1:length(cfstoplotinfofor)
                cftoplotcurr   = cfstoplotinfofor(plottingid);
                tstoplotcurr   = nsyns_tot_for_cfs_over_timesteps(find(nsyns_tot_for_cfs_over_timesteps(:,1)==cftoplotcurr),2);
                nstoplotcurr   = nsyns_tot_for_cfs_over_timesteps(find(nsyns_tot_for_cfs_over_timesteps(:,1)==cftoplotcurr),3);
                id_t_init      = find(tstoplotcurr == min(tstoplotcurr));
                id_t_final     = find(tstoplotcurr == max(tstoplotcurr));
                n_init         = nstoplotcurr(id_t_init);
                n_final        = nstoplotcurr(id_t_final);
                plot(tstoplotcurr,nstoplotcurr);
                if plottingid == 1
                    hold on;
                end
                if length(tstoplotcurr) > 1 % To make sure that all synapses
                    % made by a single CF didn't go away in one time step
                    %(which could happen if all its connections are 1 synapse
                    % big and go away at once)
                    % If there is only one time step where a climbing fiber
                    % places its synapses, then we can't fit a line to it
                    fitcurr = fit(tstoplotcurr,nstoplotcurr,'poly1');
                    if hasconverged == 1
                        avgrateofsynchangepercfperiter = [avgrateofsynchangepercfperiter; [itercurr cftoplotcurr fitcurr.p1 n_init n_final]];
                    end
                end
            end
            hold off;
            
            
        end
        
        groupedrepeatiter_premoval   = [groupedrepeatiter_premoval; p_rem_curr];
        groupedrepeatiter_powtarprob = [groupedrepeatiter_powtarprob; power_targ_curr];
        groupedrepeatiter_meanhascon = [groupedrepeatiter_meanhascon; mean(hasconvforparam)];
        groupedrepeatiter_meantcon   = [groupedrepeatiter_meantcon; mean(tconvforparam)];

        % UNCOMMENT THIS BLOCK OF CODE FOR QUANTILE ANALYSIS
% % % %         % Save quantile values for current parameters
% % % %         quantilefilename = sprintf('time-ev_and_obs_P7_quantile_results_NEW_REM_MODEL_prem_%d_powtargp_%d_niter_%d.mat',p_rem_curr,power_targ_curr,N_ITER_PER_PARAM_SET);
% % % %         quantpathandfile = ['I:\P3_P7_Comparisons\Time_Evolution_Script_for_Pub\' quantilefilename];
% % % %         save(quantpathandfile,'cumstocalculate','quantilesalliter','quantilesP7');
% % % %         
% % % %         quantile_p_vals  = []; % Columns left to right are p-values
% % % %         % for the quantiles in ascending order
% % % %         
% % % %         % Make plots that compare the quantiles for convergent time-evolved
% % % %         % connectivity matrices and the observed P7 connectivity
% % % %         for i = 1:length(cumstocalculate)
% % % %             cumcurr = cumstocalculate(i);
% % % %             cumstringaspercentage = num2str(round(cumcurr*100.0));
% % % %             valsforhistcurr       = quantilesalliter(:,i);
% % % %             edgescurr = [-0.5:1:(max(valsforhistcurr)+1.5)];
% % % %             figure()
% % % %             hcurr = histogram(valsforhistcurr,edgescurr);
% % % %             hold on;
% % % %             plot(repmat(quantilesP7(i),1,2),[0 max(hcurr.Values)],'g');
% % % %             hold off;
% % % %             title(sprintf('%s%% Quantile, P3 Time-Ev vs. P7 Obs.',cumstringaspercentage));
% % % %             xlabel('Position of Quantile');
% % % %             ylabel('Frequency');
% % % %             axis([0 max(hcurr.BinEdges)+1.0 0 max(hcurr.Values)+1.0]);
% % % %             plotnamecurrfig = sprintf('%s_pc_Quant_Time-Ev_vs_P7_Obs_prem_%d_powtargp_%d_niter_%d.fig',cumstringaspercentage,p_rem_curr,power_targ_curr,N_ITER_PER_PARAM_SET);
% % % %             plotnamecurrpng = sprintf('%s_pc_Quant_Time-Ev_vs_P7_Obs_prem_%d_powtargp_%d_niter_%d.png',cumstringaspercentage,p_rem_curr,power_targ_curr,N_ITER_PER_PARAM_SET);
% % % %             saveas(gcf,['I:\P3_P7_Comparisons\Time_Evolution_Script_for_Pub\' plotnamecurrfig]);
% % % %             saveas(gcf,['I:\P3_P7_Comparisons\Time_Evolution_Script_for_Pub\' plotnamecurrpng]);
% % % %             
% % % %             mean_for_p_val  = mean(valsforhistcurr);
% % % %             dist_obs_curr   = abs(quantilesP7(i) - mean_for_p_val);
% % % %             dists_sim_curr  = abs(valsforhistcurr - mean_for_p_val);
% % % %             p_quant_curr    = length( find(dists_sim_curr >= dist_obs_curr) )./length(valsforhistcurr);
% % % %             quantile_p_vals = [quantile_p_vals p_quant_curr];
% % % %         end
% % % %         
% % % %         quant_p_vals_fname  = sprintf('Quant_Time-Ev_vs_P7_Obs_p-vals_TS_prem_%d_powtargp_%d_niter_%d.mat',p_rem_curr,power_targ_curr,N_ITER_PER_PARAM_SET);
% % % %         quant_p_vals_file   = ['I:\P3_P7_Comparisons\Time_Evolution_Script_for_Pub\' quant_p_vals_fname];
% % % %         save(quant_p_vals_file,'cumstocalculate','quantile_p_vals');

        
    end
end

% Save the slopes and initial and final total numbers of synapses for all
% climbing fibers for all iterations in a .mat file
% save('I:\P3_P7_Comparisons\Test_Evolve_P3_to_P7_Best_Parameter_Time_and_Syn_Change_Rate_Results\test_evolve_P3_to_P7_NEW_REM_MODEL_ORIG_ADD_MDL_OPTIM_TS_prem_0p01_gam_1p3_Niter_140_avgratesofsynchangepercfperiter.mat','avgrateofsynchangepercfperiter');

% Check results
binedgesforviewing = [-0.5:1:100.5];

figure()
h_P3 = histogram(allPCconnectivitynozeros_i,binedgesforviewing);
title('P3')
xlabel('Number of Synapses');
ylabel('Frequency');
axis([0 100 0 100]);
set(gca,'FontSize',18);

figure()
h_P7 = histogram(allPCconnectivitynozeros_P7,binedgesforviewing);
title('P7')
xlabel('Number of Synapses');
ylabel('Frequency');
axis([0 100 0 100]);
set(gca,'FontSize',18);

figure()
h_f = histogram(connectivity_curr_nozeros,binedgesforviewing);
title('P3, Time-Evolved')
xlabel('Number of Synapses');
ylabel('Frequency');
axis([0 100 0 100]);
set(gca,'FontSize',18);

% Make plots in the same style used for P3 and P7 connectivity histograms
normfactP3 = sum(h_P3.Values).*h_P3.BinWidth;
normfact_f = sum(h_f.Values).*h_f.BinWidth;
normfactP7 = sum(h_P7.Values).*h_P7.BinWidth;

P3nozerosnorm     = (h_P3.Values)./normfactP3;
fnozerosnorm      = (h_f.Values)./normfact_f;
P7nozerosnorm     = (h_P7.Values)./normfactP7;

edgesP3  = h_P3.BinEdges;
edgesf   = h_f.BinEdges;
edgesP7  = h_P7.BinEdges;

bincentersP3   = 0.5 * ( edgesP3(2:end) + edgesP3(1:(end-1)) );
bincentersf    = 0.5 * ( edgesf(2:end)  + edgesf(1:(end-1))  );
bincentersP7   = 0.5 * ( edgesP7(2:end) + edgesP7(1:(end-1)) );

b1halfwidth = h_P3.BinWidth./2;
b2halfwidth = h_f.BinWidth./2;
b3halfwidth = h_P7.BinWidth./2;

figure()
b1 = createPatches(bincentersP3,P3nozerosnorm*100,b1halfwidth,'b',0.4);
axis([-0.5 70.5 0 100]);
title('P3');
xlabel('Number of Synapses Per CF-PC Pair');
ylabel('Percentage');
set(gca,'FontSize',18);

figure()
b2 = createPatches(bincentersf,fnozerosnorm*100,b2halfwidth,'b',0.4);
axis([-0.5 70.5 0 100]);
title('Simulation at P7-Like State');
xlabel('Number of Synapses Per CF-PC Pair');
ylabel('Percentage');
set(gca,'FontSize',18);


figure()
b3 = createPatches(bincentersP7,P7nozerosnorm*100,b3halfwidth,'b',0.4);
axis([-0.5 70.5 0 100]);
title('P7');
xlabel('Number of Synapses Per CF-PC Pair');
ylabel('Percentage');
set(gca,'FontSize',18);



% Inspect the average slope of the number of synapses total per climbing
% fiber as a function of time for all iterations
edgesforvisualizingavgratesofchange = [-0.0025:0.005:0.1025]; 
figure()
havg = histogram(avgrateofsynchangepercfperiter(:,3),edgesforvisualizingavgratesofchange);
title('Average Rate of Change in N Syns Per CF, All Iterations Combined');
xlabel('Avg Net Change in N Syns Per Timestep');
ylabel('Frequency');
axis([ (min(edgesforvisualizingavgratesofchange)-0.0025) (max(edgesforvisualizingavgratesofchange)+0.0025) 0 (max(havg.Values)+5) ]);


% Look at the distribution of the number of timesteps that are required for
% the simulation to converge when it converges
rowsthatconverged = find(alltrialshasconverged ==1);
ntimestepsforconvergentruns =  alltrialsnstepstoconv(rowsthatconverged);
figure()
histogram(ntimestepsforconvergentruns);
title('Number of Timesteps to Convergence for Convergent Runs');
xlabel('Number of Timesteps');
ylabel('Frequency');

% % % % % Save the times to convergence for all iterations in a .mat file
% % % % save('I:\P3_P7_Comparisons\Test_Evolve_P3_to_P7_Best_Parameter_Time_and_Syn_Change_Rate_Results\test_evolve_P3_to_P7_NEW_REM_MODEL_ORIG_ADD_MODEL_prem_0p01_powtarg_1p3_Niter_140_timestoconvergencealliters.mat','ntimestepsforconvergentruns');
% % % % 
% % % % % Save the slopes for all CFs for all iterations in a .mat file
% % % % save('I:\P3_P7_Comparisons\Test_Evolve_P3_to_P7_Best_Parameter_Time_and_Syn_Change_Rate_Results\test_evolve_P3_to_P7_NEW_REM_MODEL_ORIG_ADD_MODEL_prem_0p01_powtarg_1p3_Niter_140_avgratesofsynchangepercfperiter.mat','avgrateofsynchangepercfperiter');

return;

% Information to compute if more than one set of parameters has been
% tested

% Save results from current run into .mat files
save('I:\P3_P7_Comparisons\Test_Evolve_P3_to_P7_Parameter_Test_Results\testev_NEW_REM_MODEL_ORIG_ADD_MODEL_param_test_pr0_0p01_powtp0_1_0p1_2_0p5_4_10itpparam_hascon_tcon_all_iter_ind_next_conds_8.mat','alltrialspremoval','alltrialspowtarprob','alltrialshasconverged','alltrialsnstepstoconv');
save('I:\P3_P7_Comparisons\Test_Evolve_P3_to_P7_Parameter_Test_Results\testev_NEW_REM_MODEL_ORIG_ADD_MODEL_param_test_pr0_0p01_powtp0_1_0p1_2_0p5_4_avghascon_avgtcon_perparam_next_conds_8.mat','groupedrepeatiter_premoval','groupedrepeatiter_powtarprob','groupedrepeatiter_meanhascon','groupedrepeatiter_meantcon');

% Plot the results of the parameter search
XLIN               = linspace(min(groupedrepeatiter_premoval),max(groupedrepeatiter_premoval),30);
YLIN               = linspace(min(groupedrepeatiter_powtarprob),max(groupedrepeatiter_powtarprob),30);
[XTOPLOT,YTOPLOT]  = meshgrid(XLIN,YLIN);
MEANHASCONTOPLOT   = griddata(groupedrepeatiter_premoval,groupedrepeatiter_powtarprob,groupedrepeatiter_meanhascon,XTOPLOT,YTOPLOT,'cubic');
MEANCONVTIMETOPLOT = griddata(groupedrepeatiter_premoval,groupedrepeatiter_powtarprob,groupedrepeatiter_meantcon,XTOPLOT,YTOPLOT,'cubic');

figure()
mesh(XTOPLOT,YTOPLOT,MEANHASCONTOPLOT);
title('Mean Whether Has Convergence');
hold on;
plot3(groupedrepeatiter_premoval,groupedrepeatiter_powtarprob,groupedrepeatiter_meanhascon,'.','MarkerSize',15);
hold off;

% Re-do with filled in plot
figure()
contourf(XTOPLOT,YTOPLOT,MEANHASCONTOPLOT);
title('Mean Whether Has Convergence');
hold on;
plot3(groupedrepeatiter_premoval,groupedrepeatiter_powtarprob,groupedrepeatiter_meanhascon,'.','MarkerSize',15);
hold off;

figure()
mesh(XTOPLOT,YTOPLOT,MEANCONVTIMETOPLOT);
title('Mean Time to Convergence');
hold on;
plot3(groupedrepeatiter_premoval,groupedrepeatiter_powtarprob,groupedrepeatiter_meantcon,'.','MarkerSize',15);
hold off;

% Re-do with filled in plot
figure()
contourf(XTOPLOT,YTOPLOT,MEANCONVTIMETOPLOT);
title('Mean Time to Convergence');
hold on;
plot3(groupedrepeatiter_premoval,groupedrepeatiter_powtarprob,groupedrepeatiter_meantcon,'.','MarkerSize',15);
hold off;