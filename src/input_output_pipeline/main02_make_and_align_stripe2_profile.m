% Script to generate new, registered AP axis aligning fixed and live data
function main02_make_and_align_stripe2_profile(project,varargin)
% set default parameters
mRNA_lifetime = 7;
max_time = 40;
lambda = 1./((mRNA_lifetime)*60); % range of decay rates to calculate
ap_vec = 100:900;

for i = 1:numel(varargin)
    if ischar(varargin{i})
        eval([varargin{i} ' = varargin{i+1};'])
    end
end

% File Paths
dataPath = ['../../dat/figure_data/' project '/'];
writePath = ['../../out/input_output/' project '/'];
mkdir(writePath);
% fixed data from Dubuis 2013
fixedEvePath = '../../dat/external_data_sources/MarielaData/mov_raw_dorsal_130216_Eve_Run_Prd_orer.mat';
load(fixedEvePath);
% Load on and off time info 
load([dataPath '/inference_traces_mHMMeve2_weka_inf_2018_05_07_dT20.mat'])
load([dataPath 'analysis_data_mRNA.mat']);
%%% Obtain mapping from relative to absolute AP
ap_abs = [trace_struct_final.ap_vector];
ap_rel = [trace_struct_final.rel_ap_vector];
beta = mvregress([ones(size(ap_rel')) ap_rel'], ap_abs');
MS2_times = InterpGrid(InterpGrid<=max_time*60-70)+70;
t_vec = max_time*60-MS2_times;
%%% First Generate Predicted Profiles Corresponding to SS Model
nc_ap_array_14 = 1000*(nc_ap_array(:,nc14_vec==1)*beta(2) + beta(1));
nc_act_array_full_14 = nc_act_array_full(:,nc14_vec==1);
nc_act_array_bin_14 = nc_act_array_bin(:,nc14_vec==1);
nc_act_array_mf_14 = nc_act_array_mf(:,nc14_vec==1);


% full profile
eve2_MS2_mean_rate = NaN(numel(MS2_times),numel(ap_vec));
% binary control
eve2_MS2_mean_bin = NaN(numel(MS2_times),numel(ap_vec));
% mean rate
eve2_MS2_mean_ss = NaN(numel(MS2_times),numel(ap_vec));
% mean rate x binary
eve2_MS2_mean_ssb = NaN(numel(MS2_times),numel(ap_vec));
   
lambda_kernel = flipud(exp(-lambda*t_vec'));
% filter for nuclei that met criteria for number of frames observed during
% nc14

nc_fit_grid_full = nc_act_array_full_14(InterpGrid<=max_time*60-70,:);
nc_fit_grid_bin = nc_act_array_bin_14(InterpGrid<=max_time*60-70,:);
nc_fit_grid_mf = nc_act_array_mf_14(InterpGrid<=max_time*60-70,:);
ap_ref_array = round(nc_ap_array_14(InterpGrid<=max_time*60-70,:)*2)/2;
% make case-specific ap ref arrays
% binary
ap_ref_array_bin = ap_ref_array;
ap_ref_array_bin(isnan(nc_fit_grid_bin)) = NaN;
% full
ap_ref_array_full = ap_ref_array;
ap_ref_array_full(isnan(nc_fit_grid_full)) = NaN;
% mean rate
ap_ref_array_mf = ap_ref_array;
ap_ref_array_mf(isnan(nc_fit_grid_mf)) = NaN;

% scenario-specific ap arrays
ap_pd_rates_full = NaN(size(nc_fit_grid_full,1),numel(ap_vec));
ap_pd_rates_bin = NaN(size(nc_fit_grid_full,1),numel(ap_vec));    
ap_pd_rates_mf = NaN(size(nc_fit_grid_full,1),numel(ap_vec));    
for a = 1:numel(ap_vec)
    avg_array_full = NaN(size(nc_fit_grid_full));
    avg_array_bin = NaN(size(nc_fit_grid_full));
    avg_array_mf = NaN(size(nc_fit_grid_full));
    % random sample of indices
    % binary
    ap_indices_bin = find(ap_ref_array_bin==ap_vec(a));
    % full
    ap_indices_full = find(ap_ref_array_full==ap_vec(a));
    % mean rate
    ap_indices_mf = find(ap_ref_array_mf==ap_vec(a));    

    [row_bin, ~] = ind2sub(size(avg_array_full),ap_indices_bin);
    [row_full, ~] = ind2sub(size(avg_array_full),ap_indices_full);
    [row_mf, ~] = ind2sub(size(avg_array_full),ap_indices_mf);

    ap_pd_rates_full(:,a) = accumarray(row_full,nc_fit_grid_full(ap_indices_full),[size(nc_fit_grid_full,1) 1],@mean);
    ap_pd_rates_bin(:,a) = accumarray(row_bin,nc_fit_grid_bin(ap_indices_bin),[size(nc_fit_grid_bin,1) 1],@mean);                
    ap_pd_rates_mf(:,a) = accumarray(row_mf,nc_fit_grid_mf(ap_indices_mf),[size(nc_fit_grid_mf,1) 1],@mean);                
end
%%% apply mild spatial smoothing kernel
ap_sigma = 5; %.5% ap

ap_pd_rates_full_smooth = NaN(size(ap_pd_rates_full));
ap_pd_rates_bin_smooth = NaN(size(ap_pd_rates_full));
ap_pd_rates_mf_smooth = NaN(size(ap_pd_rates_full));

ap_mat = repmat(ap_vec,size(ap_pd_rates_full,1),1);
for i = 1:size(ap_pd_rates_full,2)
    ap_diffs = ap_mat - ap_vec(i);
    ap_weights = exp(-.5*(ap_diffs/ap_sigma).^2);
    % full
    full_vec = nansum(ap_weights.*ap_pd_rates_full,2) ./ nansum(ap_weights.*~isnan(ap_pd_rates_full),2);
    ap_pd_rates_full_smooth(:,i) = full_vec;
    % bin
    bin_vec = nansum(ap_weights.*ap_pd_rates_bin,2) ./ nansum(ap_weights.*~isnan(ap_pd_rates_bin),2);
    ap_pd_rates_bin_smooth(:,i) = bin_vec;
    % fmf
    mf_vec = nansum(ap_weights.*ap_pd_rates_mf,2) ./ nansum(ap_weights.*~isnan(ap_pd_rates_mf),2);
    ap_pd_rates_mf_smooth(:,i) = mf_vec;
end

norm_ref = conv2(ones(size(ap_pd_rates_mf_smooth)),lambda_kernel,'full');   
norm_ref_flat =  conv2(ones(size(ap_pd_rates_mf_smooth)),ones(size(lambda_kernel)),'full');   
eve_raw_mat = conv2(ap_pd_rates_full_smooth,lambda_kernel,'full')./norm_ref;
eve_raw_mat_bin = conv2(ap_pd_rates_bin_smooth,lambda_kernel,'full')./norm_ref;        
eve_raw_mat_mf = conv2(ap_pd_rates_mf_smooth,ones(size(lambda_kernel)),'full')./norm_ref_flat;        

eve2_MS2_array_full = eve_raw_mat(1:numel(lambda_kernel),:);       
eve2_MS2_array_bin = eve_raw_mat_bin(1:numel(lambda_kernel),:);
eve2_MS2_array_mf = eve_raw_mat_mf(1:numel(lambda_kernel),:);
%%% Now Calculate Offset Between Live and Fixed Gregor Data and apply
%%% correction

% find eve2 centers at 40 min 
eve2_bounds = 380:480;
[~, ind_40] = min(abs(MS2_times-40*60));
model_profile = eve2_MS2_array_full(ind_40,:);
fixed_profile = EveMov_final(40,eve2_bounds);
[~, mi_model] = nanmax(model_profile);
mi_model = mi_model + min(ap_vec);
[~, mi_fixed] = nanmax(fixed_profile);
mi_fixed = mi_fixed + min(eve2_bounds);
% calculate correction
ap_shift = mi_fixed - mi_model;

%%% save
%%% time-dependent profiles
eve2_struct.eve2_MS2_mean_full = eve2_MS2_array_full;
eve2_struct.eve2_MS2_mean_norm_full = (eve2_MS2_array_full-min(eve2_MS2_array_full(:)))./(max(eve2_MS2_array_full(:))-min(eve2_MS2_array_full(:)));

eve2_struct.eve2_MS2_mean_bin = eve2_MS2_array_bin;
eve2_struct.eve2_MS2_mean_norm_bin = (eve2_MS2_array_bin-min(eve2_MS2_array_bin(:)))./(max(eve2_MS2_array_bin(:))-min(eve2_MS2_array_bin(:)));

eve2_struct.eve2_MS2_mean_mf = eve2_MS2_array_mf;
eve2_struct.eve2_MS2_mean_norm_mf = (eve2_MS2_array_mf-min(eve2_MS2_array_mf(:)))./(max(eve2_MS2_array_mf(:))-min(eve2_MS2_array_mf(:)));

eve2_struct.lambda = lambda;
eve2_struct.MS2_times = MS2_times;
eve2_struct.ap_shift = ap_shift;
save([writePath 'eve2_struct.mat'],'eve2_struct')
