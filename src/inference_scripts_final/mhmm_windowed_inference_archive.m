% Script to Conduct Windowed HMM Inference on Experimental Data
function mhmm_windowed_inference_archive(project,outPath,inference_data,bin_groups,id_string,...
    Tres,K,w,t_window,inference_times,varargin)

% set defaults
% inference parameters
n_localEM = 25; % set num local runs
n_steps_max = 500; % set max steps per inference
eps = 1e-4; % set convergence criteria
bootstrap = 1;
n_bootstrap = 10;
sample_size = 8000;
min_dp_per_inf = 1250; % inference will be aborted if fewer present   
% assorted qc variables
min_dp = 10; % min length of traces to include
clipped = 1; % if 0 use "full" trace with leading and trailing 0's
clipped_ends = 1; % if one, remove final w time steps from traces
                
savio = 0; % set to 1 if inference is being conducted on savio cluster

for i = 1:numel(varargin)
    if ischar(varargin{i})
        eval([varargin{i} ' = varargin{i+1};'])
    end
end

addpath('../utilities'); % Route to utilities folder
warning('off','all') %Shut off Warnings                                
% max num workers
if savio
    MaxWorkers = 24;
else
    MaxWorkers = 25;
end
if bootstrap
    d_type = '_boot';
end
alpha = inference_data(1).alpha_frac*w; % Rise Time for MS2 Loops
% Set write path (inference results are now written to external directory)
out_suffix =  ['/' project '/' id_string '/']; 
out_dir = [outPath out_suffix];
mkdir(out_dir);

%%%Set write directories
if clipped_ends
    end_clip = w + 1;
else
    end_clip = 0;
end
% apply time filtering 
trace_struct_filtered = [];
for i = 1:length(inference_data)
    temp = inference_data(i);   
    time = temp.time_interp(w+1:end-end_clip)/60; % we must ignore first w + 1 time points for windowed inference
    fluo = temp.fluo_interp(w+1:end-end_clip);               
    
    if numel(time) >= min_dp
        temp.fluo = fluo;
        temp.time = time;
        trace_struct_filtered = [trace_struct_filtered temp];
    end
end

trace_struct_filtered = trace_struct_filtered([trace_struct_filtered.inference_flag]==1);
% get last obs time for each set
set_vec_trace = [trace_struct_filtered.setID];
set_index = unique(set_vec_trace);

% trace start and end times for each trace
first_time_vec = [];
last_time_vec = [];
for i = 1:length(trace_struct_filtered)
    first_time_vec = [first_time_vec min(trace_struct_filtered(i).time)];
    last_time_vec = [last_time_vec max(trace_struct_filtered(i).time)];
end
%% Conduct Inference
% structure array to store the analysis data
local_meta = struct; % Store local_EM results

for g = 1:length(bin_groups) % loop through different AP groups
    bin_list = bin_groups{g}; % get groups for present iteration         
    for t = 1:length(inference_times)
        t_inf = inference_times(t);
        t_start = t_inf - t_window/2;
        t_stop = t_inf + t_window/2;              
        for b = 1:n_bootstrap
            iter_start = now;
            s = (g-1)*n_bootstrap + b;
            local_struct = struct;
            init_struct = struct;
            output = struct;           
            % Use current time as unique inference identifier 
            inference_id = num2str(round(10e5*now));
            % Generate filenames           
            fName_sub = ['eveSet_w' num2str(w) '_bin' num2str(round(10*mean(bin_list))) '_K' num2str(K) '_t' inference_id];                
            out_file = [out_dir '/' fName_sub];
            % Sample with replacement from available sets
            set_boot_vec = randsample(set_index,numel(set_index),true);
            % Now sample trace_struct
            trace_struct_boot = [];
            first_time_boot = [];
            last_time_boot = [];
            for i = 1:numel(set_boot_vec)
                trace_struct_boot = [trace_struct_boot trace_struct_filtered(set_vec_trace==set_boot_vec(i))];
                first_time_boot = [first_time_boot first_time_vec(set_vec_trace==set_boot_vec(i))];
                last_time_boot = [last_time_boot last_time_vec(set_vec_trace==set_boot_vec(i))];
            end                       
            % generate time-resolved reference vector
            ap_ref_vec = NaN(1,length(trace_struct_boot));
            for i = 1:length(trace_struct_boot)
                tt = trace_struct_boot(i).time;
                ap = trace_struct_boot(i).rel_ap_vector_interp;            
                tr_ap = round(mean(ap(tt>=t_start&tt<t_stop)));            
                ap_ref_vec(i) = tr_ap;
            end  
            % Extract fluo_data
            trace_filter = ismember(ap_ref_vec,bin_list) & (last_time_boot-min_dp*Tres/60) >= t_start ...
                    & (first_time_boot+min_dp*Tres/60) < t_stop; 
            trace_ind = find(trace_filter);
            % make final inference set
            inference_set = [];
            for i = trace_ind
                temp = trace_struct_boot(i);
                tt = temp.time;
                ft = temp.fluo;
                temp.time = tt(tt>=t_start & tt < t_stop);
                temp.fluo = ft(tt>=t_start & tt < t_stop);
                if numel(temp.fluo) > min_dp && sum(temp.fluo>0) > 1 % exclude strings of pure zeros
                    inference_set = [inference_set temp];                    
                end
            end
            if isempty(inference_set)
                set_size = 0;
                continue
            else
                set_size = length([inference_set.fluo]);                 
            end
            skip_flag = set_size < min_dp_per_inf ;
            sample_index = 1:length(inference_set);
            if bootstrap                        
                ndp = 0;    
                sample_ids = [];                    
                %Reset bootstrap size to be on order of set size for small bins
                sample_size = min(sample_size,ceil(set_size/100)*100);                    
                while ndp < sample_size
                    tr_id = randsample(sample_index,1);
                    sample_ids = [sample_ids tr_id];
                    ndp = ndp + length(inference_set(tr_id).time);
                end
                fluo_data = cell([length(sample_ids), 1]);    
                time_data = cell([length(sample_ids), 1]);    
                sample_particles = [inference_set(sample_ids).ParticleID];
                for tr = 1:length(sample_ids)
                    fluo_data{tr} = inference_set(sample_ids(tr)).fluo;                    
                    time_data{tr} = inference_set(sample_ids(tr)).time;                    
                end            
            else % Take all relevant traces if not bootstrapping
                fluo_data = cell([length(inference_set), 1]);            
                for tr = 1:length(inference_set)
                    fluo_data{tr} = inference_set(tr).fluo;
                    time_data{tr} = inference_set(tr).time;                    
                end
            end
            % Random initialization of model parameters
            param_init = initialize_random (K, w, fluo_data);
            % Approximate inference assuming iid data for param initialization                
            local_iid_out = local_em_iid_reduced_memory_truncated (fluo_data, param_init.v, ...
                                param_init.noise, K, w, alpha, n_steps_max, eps);
            noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
            v_iid = exp(local_iid_out.v_logs);            
            p = gcp('nocreate');
            if isempty(p)
                parpool(MaxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
            elseif p.NumWorkers > MaxWorkers
                delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                parpool(MaxWorkers);
            end
            parfor i_local = 1:n_localEM % Parallel Local EM                
                % Random initialization of model parameters
                param_init = initialize_random_with_priors(K, noise_iid, v_iid);
                % Get Intial Values
                pi0_log_init = log(param_init.pi0);
                A_log_init = log(param_init.A);
                v_init = param_init.v;                        
                noise_init = param_init.noise;
                % Record
                init_struct(i_local).A_init = exp(A_log_init);                
                init_struct(i_local).v_init = v_init;
                init_struct(i_local).noise_init = noise_init;                
                init_struct(i_local).subset_id = i_local;
                %--------------------LocalEM Call-------------------------%
                local_out = local_em_MS2_reduced_memory_truncated(fluo_data, ...
                    v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
                    alpha, n_steps_max, eps);                    
                %---------------------------------------------------------%                
                % Save Results 
                local_struct(i_local).inference_id = inference_id;
                local_struct(i_local).subset_id = i_local;
                local_struct(i_local).logL = local_out.logL;                
                local_struct(i_local).A = exp(local_out.A_log);
                local_struct(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
                local_struct(i_local).r = exp(local_out.v_logs).*local_out.v_signs / Tres;                                
                local_struct(i_local).noise = 1/exp(local_out.lambda_log);
                local_struct(i_local).pi0 = exp(local_out.pi0_log);
%                         local_struct(i_local).total_time = local_out.runtime;
                local_struct(i_local).total_steps = local_out.n_iter;               
                local_struct(i_local).soft_struct = local_out.soft_struct;               
            end
            local_meta(s).init = init_struct;
            local_meta(s).local = local_struct;
            [logL, max_index] = max([local_struct.logL]); % Get index of best result                    
            % Save parameters from most likely local run
            output.pi0 =local_struct(max_index).pi0;                        
            output.r = local_struct(max_index).r(:);           
            output.noise = local_struct(max_index).noise;
            output.A = local_struct(max_index).A(:);
            output.A_mat = local_struct(max_index).A;            
            % get soft-decoded structure
            output.soft_struct = local_struct(max_index).soft_struct;
            % Info about run time
            output.total_steps = local_struct(max_index).total_steps;                                  
            output.total_time = 100000*(now - iter_start);            
            % Save inference ID variables
            output.APbin = min(bin_list):max(bin_list);
            % other inference characteristics
            output.t_window = t_window;
            output.t_inf = t_inf;                
            output.bootstrap_flag = bootstrap;   
            output.iter_id = b;
            output.start_time_inf = 0;                              
            output.clipped = clipped;
            output.clipped_ends = clipped_ends;
            if bootstrap || set_bootstrap                    
                output.particle_ids = sample_particles;                
                output.set_boot_vec = set_boot_vec;
            end                
            output.N = ndp;
            output.w = w;
            output.K = K;
            output.alpha = alpha;
            output.deltaT = Tres; 
            % save inference data used
            output.fluo_data = fluo_data;
            output.time_data = time_data;
            % record whether current inference has too few data points to
            % use
            output.skip_flag = skip_flag;
            save([out_file '.mat'], 'output');           
        end  
    end
end    
