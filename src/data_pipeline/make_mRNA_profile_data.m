% Script to Generate mRNA Profiles to Represent Different Models Discussed
% in Paper

function main03_make_mRNA_profile_data(project,varargin)

% indicate whether to use nc14 QC vec
use_nc14_qc_vec = 1;
n_boots = 100; % number of bootstraps to use for error estimation
max_time = 50; % desired ref time (in minutes)
res = .5; % AP resolution
mRNA_lifetime_vec = [1 7:4:15 21 Inf]; % mRNA lifetimes to simulate
for i = 1:numel(varargin)
    if ischar(varargin{i})
        if ismember(varargin{i},{'n_boots','res','max_time','use_nc14_qc_vec'})
            eval([varargin{i} ' = varargin{i+1};'])
        end
    end
end
ap_plot = -8:res:8;
DataPath = ['../../dat/figure_data/' project '/'];

% Load on and off time info 
load([DataPath 'analysis_data_mRNA.mat']);
DataPath = ['../../dat/figure_data/' project '/'];
% if you do not want to use nc14 qc, overwrite the vector
save_suffix = 'qc1';    
if ~use_nc14_qc_vec
    save_suffix = 'qc0';
    nc14_vec = true(size(nc14_vec));
end

if strcmp(project,'eve7stripes_inf_2018_04_28')
    save_suffix = [save_suffix '_BAC'];
elseif strcmp(project,'mHMMeve2_weka_inf_2018_05_07')
    save_suffix = [save_suffix '_reporter'];
end

% filter arrays
nc_set_vec = NucSets(nc14_vec==1);
set_index = unique(nc_set_vec);
nc_ap_array_14 = nc_ap_array(:,nc14_vec==1);
nc_act_array_full_14 = nc_act_array_full(:,nc14_vec==1);
nc_act_array_bin_14 = nc_act_array_bin(:,nc14_vec==1);
nc_act_array_mf_14 = nc_act_array_mf(:,nc14_vec==1);
% filter for nuclei that met criteria for number of frames observed during
 % nc14
nc_fit_grid_full = nc_act_array_full_14(InterpGrid<=max_time*60-70,:);
nc_fit_grid_bin = nc_act_array_bin_14(InterpGrid<=max_time*60-70,:);
nc_fit_grid_mf = nc_act_array_mf_14(InterpGrid<=max_time*60-70,:);
ap_ref_array = round(nc_ap_array_14(InterpGrid<=max_time*60-70,:)*1/res)*res;

% generate on and off time vectors
on_time_vec = NaN(1,size(nc_fit_grid_bin,2));
off_time_vec = NaN(1,size(nc_fit_grid_bin,2));
on_ap_vec = NaN(1,size(nc_fit_grid_bin,2));
off_ap_vec = NaN(1,size(nc_fit_grid_bin,2));
for i = 1:size(nc_fit_grid_bin,2)
    act_indices = find(nc_fit_grid_bin(:,i)>0);
    if ~isempty(act_indices)
        on_time_vec(i) = InterpGrid(act_indices(1));
        on_ap_vec(i) = ap_ref_array(act_indices(1),i);
        off_time_vec(i) = InterpGrid(act_indices(end));
        off_ap_vec(i) = ap_ref_array(act_indices(end),i);
    end
end

% ap reference array
ap_ref_array_full = ap_ref_array;
ap_ref_array_full(isnan(nc_fit_grid_full)) = NaN;
% indicates nuclei that turned on at some point
ever_on_vec = nanmax(nc_fit_grid_full)>0;
duration_vec = nansum(nc_fit_grid_bin);
lifetime_vec = sum(~isnan(nc_fit_grid_full));
% arrays to store metrics
fraction_on_array = NaN(n_boots,length(ap_plot)); % fraction ever active
ss_mRNA_rate_array = NaN(n_boots,length(ap_plot)); % effect of switching alon
duration_array = NaN(n_boots,length(ap_plot)); % effective duration of activity
on_time_array = NaN(n_boots,length(ap_plot)); % average on times
off_time_array = NaN(n_boots,length(ap_plot)); % average on times


for n = 1:n_boots
    % bootstrap sets first 
    boot_sets = randsample(set_index,numel(set_index),true);
    available_ids = [];
    for s = 1:numel(boot_sets)
        available_ids = [available_ids find(nc_set_vec==boot_sets(s))];
    end
    boot_nc_ids = randsample(available_ids,numel(nc_set_vec),true);
    ap_ref_boot = ap_ref_array_full(:,boot_nc_ids);
    mf_grid_boot = nc_fit_grid_mf(:,boot_nc_ids);
    ever_on_vec_boot = ever_on_vec(boot_nc_ids);
    duration_vec_boot = duration_vec(boot_nc_ids);
    lifetime_vec_boot = lifetime_vec(boot_nc_ids);
    on_time_boot = on_time_vec(boot_nc_ids);
    off_time_boot = off_time_vec(boot_nc_ids);
    on_ap_boot = round(on_ap_vec(boot_nc_ids),2);
    off_ap_boot = round(off_ap_vec(boot_nc_ids),2);
    for a = 1:length(ap_plot)
        ap = ap_plot(a);    
        ap_flag_array = round(ap_ref_boot,2)==round(ap,2);                     
        ss_mRNA_rate_array(n,a) = nanmean(mf_grid_boot(ap_flag_array));        
        
        % duration calculation. 
        ap_ft = nanmax(ap_flag_array)==1; % ever present in AP region
        ever_on = ever_on_vec_boot(ap_ft);        
        act_time_vec = duration_vec_boot(ap_ft);
        ap_time = sum(ap_flag_array(:,ap_ft));          
        total_time_vec = lifetime_vec_boot(ap_ft);        
        ap_frac_vec = ap_time ./ total_time_vec; % fraction of total lifetime in AP region
        
        on_ft = act_time_vec>0;
        duration_array(n,a) = nansum(act_time_vec(on_ft).*ap_frac_vec(on_ft))/nansum(ap_frac_vec(on_ft));        
        fraction_on_array(n,a) = nansum(ever_on.*ap_frac_vec)/nansum(ap_frac_vec);
        
        on_time_array(n,a) = nanmean(on_time_boot(on_ap_boot==round(ap,2)));
        off_time_array(n,a) = nanmean(off_time_boot(off_ap_boot==round(ap,2)));
    end   
end

duration_mean = nanmean(duration_array,1);
duration_ste = nanstd(duration_array,1);

fraction_on_mean = nanmean(fraction_on_array,1);
fraction_on_ste = nanstd(fraction_on_array,1);

ss_mRNA_mean = nanmean(ss_mRNA_rate_array,1);
ss_mRNA_ste = nanstd(ss_mRNA_rate_array,1);

ss_full_mean = nanmean(ss_mRNA_rate_array.*fraction_on_array);
ss_full_ste = nanstd(ss_mRNA_rate_array.*fraction_on_array);

ap_on_time_mean = nanmean(on_time_array);
ap_on_time_ste = nanstd(on_time_array);

ap_off_time_mean = nanmean(off_time_array);
ap_off_time_ste = nanstd(off_time_array);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now Calculate full and binary profiles--these quantities may vary over
%%% course of nc14
lambda_vec = 1./((mRNA_lifetime_vec)*60); % range of decay rates to calculate
MS2_times = InterpGrid(InterpGrid<=max_time*60-70)+70;
t_vec = max_time*60-MS2_times;
% full profile
eve2_MS2_mean_full_lambda = NaN(numel(MS2_times),numel(ap_plot),numel(lambda_vec));
eve2_MS2_ste_full_lambda = NaN(numel(MS2_times),numel(ap_plot),numel(lambda_vec));
% binary control
eve2_MS2_mean_bin_lambda = NaN(numel(MS2_times),numel(ap_plot),numel(lambda_vec));
eve2_MS2_ste_bin_lambda = NaN(numel(MS2_times),numel(ap_plot),numel(lambda_vec));
% mean rate
eve2_MS2_mean_mf_lambda = NaN(numel(MS2_times),numel(ap_plot),numel(lambda_vec));
eve2_MS2_ste_mf_lambda = NaN(numel(MS2_times),numel(ap_plot),numel(lambda_vec));
% mean rate x binary
eve2_MS2_mean_ssb_lambda = NaN(numel(MS2_times),numel(ap_plot),numel(lambda_vec));
eve2_MS2_ste_ssb_lambda = NaN(numel(MS2_times),numel(ap_plot),numel(lambda_vec));

for i = 1:numel(lambda_vec)    
    lambda = lambda_vec(i);
    lambda_kernel = flipud(exp(-lambda*t_vec'));       
    %%% Calculate average production rate per ap bin as a function of time
    eve2_MS2_array_full = NaN(size(nc_fit_grid_full,1),numel(ap_plot),n_boots);
    eve2_MS2_array_bin = NaN(size(nc_fit_grid_full,1),numel(ap_plot),n_boots);
    eve2_MS2_array_mf = NaN(size(nc_fit_grid_full,1),numel(ap_plot),n_boots);
    tic
    parfor n = 1:n_boots
        % bootstrap sets first 
        boot_sets = randsample(set_index,numel(set_index),true);
        available_ids = [];
        for s = 1:numel(boot_sets)
            available_ids = [available_ids find(nc_set_vec==boot_sets(s))];
        end
        % Take bootstrap sample of nuclei
        boot_nc_ids = randsample(available_ids,numel(nc_set_vec),true);
        ap_ref_boot = ap_ref_array_full(:,boot_nc_ids);
        mf_grid_boot = nc_fit_grid_mf(:,boot_nc_ids);
        full_grid_boot = nc_fit_grid_full(:,boot_nc_ids);
        bin_grid_boot = nc_fit_grid_bin(:,boot_nc_ids);
        % Arrays to store average "production rate" as defined by each
        % metric
        ap_pd_rates_full = NaN(size(nc_fit_grid_full,1),numel(ap_plot));
        ap_pd_rates_bin = NaN(size(nc_fit_grid_full,1),numel(ap_plot));    
        ap_pd_rates_mf = NaN(size(nc_fit_grid_full,1),numel(ap_plot));    
        % Calculate average production rate for each AP-time unit
        for a = 1:numel(ap_plot)
            ap_flag_array = double(round(ap_ref_boot,2) == round(ap_plot(a),2));            
            ap_flag_array(ap_flag_array==0) = NaN;
            ap_pd_rates_full(:,a) = nanmean(full_grid_boot.*ap_flag_array,2);
            ap_pd_rates_bin(:,a) = nanmean(bin_grid_boot.*ap_flag_array,2);
            ap_pd_rates_mf(:,a) = nanmean(mf_grid_boot.*ap_flag_array,2);
            
            ap_pd_rates_full(isnan(ap_pd_rates_full)) = 0;
            ap_pd_rates_bin(isnan(ap_pd_rates_bin)) = 0;
            ap_pd_rates_mf(isnan(ap_pd_rates_mf)) = 0;
        end
        % convolve with decay and averaging kernels as appropriate
        norm_ref_lambda = conv2(ones(size(ap_pd_rates_mf)),lambda_kernel,'full');
        norm_ref_avg = conv2(ones(size(ap_pd_rates_mf)),ones(size(lambda_kernel)),'full');
        % use convolutions to calculate temporal profiles
        eve_raw_mat_full = conv2(ap_pd_rates_full,lambda_kernel,'full')./norm_ref_lambda;
        eve_raw_mat_bin = conv2(ap_pd_rates_bin,lambda_kernel,'full')./norm_ref_lambda;        
        eve_raw_mat_mf = conv2(ap_pd_rates_mf,ones(size(lambda_kernel)),'full')./norm_ref_avg;        

        eve2_MS2_array_full(:,:,n) = eve_raw_mat_full(1:numel(lambda_kernel),:);       
        eve2_MS2_array_bin(:,:,n) = eve_raw_mat_bin(1:numel(lambda_kernel),:);
        eve2_MS2_array_mf(:,:,n) = eve_raw_mat_mf(1:numel(lambda_kernel),:);                
    end
    % Full mode
    eve2_MS2_mean = nanmean(eve2_MS2_array_full,3);
    eve2_MS2_ste = nanstd(eve2_MS2_array_full,[],3);

    eve2_MS2_mean_full_lambda(:,:,i) = eve2_MS2_mean;
    eve2_MS2_ste_full_lambda(:,:,i) = eve2_MS2_ste;
    % Binary model
    eve2_MS2_mean_bin = nanmean(eve2_MS2_array_bin,3);
    eve2_MS2_ste_bin = nanstd(eve2_MS2_array_bin,[],3);
    
    eve2_MS2_mean_bin_lambda(:,:,i) = eve2_MS2_mean_bin;
    eve2_MS2_ste_bin_lambda(:,:,i) = eve2_MS2_ste_bin;
    % MF Model
    eve2_MS2_mean_mf = nanmean(eve2_MS2_array_mf,3);
    eve2_MS2_ste_mf = nanstd(eve2_MS2_array_mf,[],3);

    eve2_MS2_mean_mf_lambda(:,:,i) = eve2_MS2_mean_mf;
    eve2_MS2_ste_mf_lambda(:,:,i) = eve2_MS2_ste_mf;
    
    % MF Model x Binary Model
    eve2_MS2_mean_ssb = nanmean(eve2_MS2_array_mf.*eve2_MS2_array_bin,3);
    eve2_MS2_ste_ssb = nanstd(eve2_MS2_array_mf.*eve2_MS2_array_bin,[],3);
    
    eve2_MS2_mean_ssb_lambda(:,:,i) = eve2_MS2_mean_ssb;
    eve2_MS2_ste_ssb_lambda(:,:,i) = eve2_MS2_ste_ssb;
end

%%% save
%%% time-dependent profiles
out_struct.eve2_MS2_mean_ssb = eve2_MS2_mean_ssb_lambda;
out_struct.eve2_MS2_ste_ssb = eve2_MS2_ste_ssb_lambda;
out_struct.eve2_MS2_mean_bin = eve2_MS2_mean_bin_lambda;
out_struct.eve2_MS2_ste_bin = eve2_MS2_ste_bin_lambda;
out_struct.eve2_MS2_mean = eve2_MS2_mean_full_lambda;
out_struct.eve2_MS2_ste = eve2_MS2_ste_full_lambda;
out_struct.eve2_MS2_mean_mf = eve2_MS2_mean_mf_lambda;
out_struct.eve2_MS2_ste_mf = eve2_MS2_ste_mf_lambda;
out_struct.lambda_vec = lambda_vec;
out_struct.MS2_times = MS2_times;

%%% SS profiles
out_struct.ss_full_mean = ss_full_mean;
out_struct.ss_full_ste = ss_full_ste;
out_struct.ss_MS2_mean = ss_mRNA_mean;
out_struct.ss_MS2_ste = ss_mRNA_ste;
out_struct.fraction_on_mean = fraction_on_mean;
out_struct.fraction_on_ste = fraction_on_ste;
out_struct.duration_mean = duration_mean;
out_struct.duration_ste = duration_ste;
out_struct.ap_index = ap_plot;

% On/Off times
out_struct.ap_on_time_mean = ap_on_time_mean;
out_struct.ap_on_time_ste = ap_on_time_ste;

out_struct.ap_off_time_mean = ap_off_time_mean;
out_struct.ap_off_time_ste = ap_off_time_ste;

save([DataPath 'profile_data_' save_suffix '.mat'],'out_struct')
