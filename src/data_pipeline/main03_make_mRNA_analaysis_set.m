% Script to Generate Stripe Profile Arrays and Extract Time On/Off info

% Assumes that Instantaneous fluorescence is a reasonable proxy for
% instantaneous productivity (lagged by w/2). See Appendix 5

% NL: Revised on 11.07.2018 for second eLife Submission

function main03_make_mRNA_analaysis_set(project,varargin)
n_boots = 100; 
ap_index = (30:50);
dataPath = '../../dat/';
InterpGrid = 0:20:60*60;
for i = 1:numel(varargin)
    if ischar(varargin{i})
        if ismember(varargin{i},{'n_boots','ap_index'})
            eval([varargin{i} ' = varargin{i+1};'])
        end
    end
end

load([dataPath '/' project '/qc_nucleus_struct.mat'])

% Generate indexing vectors and initialo
% InterpGrid = nucleus_struct(1).InterpGrid; % Time grid used for interpolation
NucSets = [nucleus_struct.setID];
% generate qc flag that filters for nuclei present for majority of nuclear
% cycle
nc14_vec = [nucleus_struct.nc_trace_on_flag];
gtype_id_vec = [];
for i = 1:numel(nucleus_struct)
    gtype_id_vec = [gtype_id_vec nucleus_struct(i).gtypeID(1)];    
end
gtype_str_vec = {nucleus_struct.genotype};
gtype_id_index = unique(gtype_id_vec,'stable');
gtype_str_index = unique(gtype_str_vec,'stable');
nGenotypes = numel(gtype_id_index);
%%% Make arrays that aggregate nuclear and spot tracking info
nc_act_array_mf = NaN(numel(InterpGrid),numel(nc14_vec)); % mean field (avg among active)
nc_act_array_full = NaN(numel(InterpGrid),numel(nc14_vec)); % full (avg among all extant nuclei)
nc_act_array_bin = NaN(numel(InterpGrid),numel(nc14_vec)); % full (avg among all extant nuclei)
nc_ap_array = NaN(numel(InterpGrid),numel(nc14_vec)); % full (avg among all extant nuclei)
% Iterate through nucleus structure
for i = 1:length(nucleus_struct)
    nc_time = nucleus_struct(i).time_interp;    
    nc_ap = 100*nucleus_struct(i).ap_vector_interp;           
    % add AP positions
    FOV_ft = ismember(InterpGrid,round(nc_time));
    nc_ap_array(FOV_ft,i) = nc_ap;  
    % in full case count nucleus in mean wherever present, regardless of whether it has a spot   
    nc_act_array_full(FOV_ft,i) = 0;  
    last_index_nc = find(FOV_ft==1,1,'last');
    % if nucleus had an active locus, add transcription info
    pt_fluo = nucleus_struct(i).fluo_interp;
    nc_act_array_mf(FOV_ft,i) = pt_fluo;
   
    f_ft = ismember(InterpGrid,nc_time(~isnan(pt_fluo)));

    nc_act_array_full(f_ft,i) = pt_fluo(~isnan(pt_fluo));
    % now update binary array
    last_index_tr = find(f_ft,1,'last');
    if ~isempty(last_index_tr)
        nc_act_array_bin(f_ft,i) = 1; % 1 while locus active       
        nc_act_array_bin(last_index_tr+1:last_index_nc,i) = 0; % 0 thereafter until nucleus leaves FOV
    end    
end

%%% Calculate average on and off times for each AP region--use average nucleus AP
on_boot_mat = NaN(length(ap_index),nGenotypes,n_boots);
off_boot_mat = NaN(length(ap_index),nGenotypes,n_boots);
ss_mean_boot_mat = NaN(length(ap_index),nGenotypes,n_boots);
for g = 1:nGenotypes
    for a = 1:length(ap_index)
        ap = ap_index(a);
        % grab relevant traces
        tr_filter = nc14_vec & round(nanmean(nc_ap_array)) == ap & gtype_id_vec == g;
        if ~any(tr_filter)
            continue
        end
        nc_act_filt = nc_act_array_mf(:,tr_filter);
        ap_traces = nc_act_array_mf(:,tr_filter);
        ap_sets = NucSets(tr_filter);
        ap_on_times = [];
        ap_off_times = [];
        for p = 1:size(ap_traces,2)
            ap_on_times = [ap_on_times InterpGrid(find(~isnan(ap_traces(:,p)),1))];
            ap_off_times = [ap_off_times InterpGrid(find(~isnan(ap_traces(:,p)),1,'last'))];
        end    

        % take bootstrap samples
        for n = 1:n_boots
            % first bootstrap at level of sets
            if numel(unique(ap_sets)) > 1
                boot_sets = randsample(unique(ap_sets),numel(unique(ap_sets)),true);
            else
                boot_sets = ap_sets;
            end
            pt_ids = [];
            for s = 1:numel(boot_sets)
                pt_ids = [pt_ids find(ap_sets==boot_sets(s))];
            end      
            boot_samp = randsample(pt_ids,length(pt_ids),true);       
            on_boot_mat(a,g,n) = mean(ap_on_times(boot_samp));
            off_boot_mat(a,g,n) = mean(ap_off_times(boot_samp));
            nc_samp = nc_act_filt(:,boot_samp);
            ss_mean_boot_mat(a,g,n) = nanmean(nc_samp(:));
        end   
    end
end

% save relevant fields analysis_data_mRNA
analysis_data_mRNA = struct;
% ap trends
analysis_data_mRNA.mean_ap_on = nanmean(on_boot_mat,3);
analysis_data_mRNA.ste_ap_on = nanstd(on_boot_mat,[],3);
analysis_data_mRNA.mean_ap_off = nanmean(off_boot_mat,3);
analysis_data_mRNA.ste_ap_off = nanstd(off_boot_mat,[],3);
analysis_data_mRNA.ss_mean_rate = nanmean(ss_mean_boot_mat,3);
analysis_data_mRNA.ss_mean_ste = nanstd(ss_mean_boot_mat,[],3);
% qc vectors
analysis_data_mRNA.nc14_vec = nc14_vec;
analysis_data_mRNA.nc_qc_flag = [nucleus_struct.nc_qc_flag];
% indexing vecotrs
analysis_data_mRNA.gtype_id_vec = gtype_id_vec;
analysis_data_mRNA.gtype_id_index = gtype_id_index;
analysis_data_mRNA.gtype_str_index = gtype_str_index;
analysis_data_mRNA.setID_vec = [nucleus_struct.setID];
% nucleus activity arrays
analysis_data_mRNA.nc_ap_array = nc_ap_array;
analysis_data_mRNA.nc_act_array_full = nc_act_array_full;
analysis_data_mRNA.nc_act_array_mf = nc_act_array_mf;
analysis_data_mRNA.nc_act_array_bin = nc_act_array_bin;
% clear unwanted variables

save([dataPath '/' project  '/analysis_data_mRNA.mat'],'analysis_data_mRNA')