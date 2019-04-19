% Script to generate regression weights to ensure that all spatiotemporal
% regions are accorded equal priority in input/output regressions
function main04_make_regression_weights(project,varargin)

% declare repsonse variables of interest 
response_var_flags = {'bin_act_flag','z_on_flag'};
response_vars = {'on','z_on'};
% set bounds and level of granularity for weight calculation
ap_index_cell = {320:10:480 320:10:480}; 
time_index_cell = {10*60:60:40*60 10*60:60:40*60};
% smoothing kernel parameters
kernel_size = 5;
kernel_sigma = 2; % 2% AP
% min data points per region to use
min_dp = 20; % cap region weights
% data paths
inputOutputPath = ['../../out/input_output/' project '/'];

for i = 1:numel(varargin)
    if ischar(varargin{i})
        eval([varargin{i} ' = varargin{i+1};'])
    end
end


% load data
regression_set = readtable([inputOutputPath 'regression_set.csv']);
% indexing vectors
ap_vec = 10*round(regression_set.APOrig/10);
time_vec = 60*round(regression_set.Time/60);

% Calculate sample weights
for f = 1:numel(response_var_flags)
    ap_index = ap_index_cell{f};
    time_index = time_index_cell{f};    
    dp_var_flag = response_var_flags{f};    
    % update flag var to account for space and time bounds 
    flag_vec = regression_set.(dp_var_flag);       
    flag_vec_new = flag_vec==1&ismember(ap_vec,ap_index)&ismember(time_vec,time_index);
    regression_set.(dp_var_flag) = flag_vec_new;

    filt_set = regression_set(regression_set.(dp_var_flag)==1,:);    
    filt_set.APOrig = round(filt_set.APOrig/10)*10;
    filt_set.Time = round(filt_set.Time/60)*60;
    % generate summary structure
    statarray = grpstats(filt_set,{'APOrig','Time'},{'mean'},'DataVars',response_vars{f});           
    % calculate reweighting factors
    n_vec = statarray.GroupCount;            
    % get raw spatio-temporal mean grid
    ct_grid = NaN(numel(time_index),numel(ap_index));
    mf_vec = statarray.(['mean_' response_vars{f}]);    
    for a = 1:numel(ap_index)
        for t = 1:numel(time_index)        
            mf = mf_vec(statarray.APOrig==ap_index(a)&statarray.Time==time_index(t));
            if ~isempty(mf)          
                ct_grid(t,a) = n_vec(statarray.APOrig==ap_index(a)&statarray.Time==time_index(t));
            end
        end
    end         
    % now use gaussian kernel to obtain smoothed, average profile while
    % maintaining initial resolution
    norm_array = NaN(size(ct_grid));
    norm_array(~isnan(ct_grid)) = 1;
    imageFilter = fspecial('gaussian',kernel_size,kernel_sigma);
    ctFiltered = nanconv(ct_grid,imageFilter, 'nonanout');      
    
    % dataFiltered = imfilter(ct_grid,imageFilter);  
    refFiltered = nanconv(norm_array,imageFilter, 'nonanout');  
    ctFilteredNorm = ctFiltered ./ refFiltered;            
    ctFilteredNorm(ctFilteredNorm<min_dp) = NaN;
    
    % make re-weighting array
    weight_array = ctFilteredNorm;  
    weight_array = max(weight_array(:)) ./ weight_array; 
    % map back to full set
    ap_vec_full = round(regression_set.APOrig/10)*10;
    time_vec_full = round(regression_set.Time/60)*60;
    wt_vec_full = NaN(size(regression_set,1),1);
    for a = 1:numel(ap_index)
        for t = 1:numel(time_index)
            wt_vec_full(ap_vec_full==ap_index(a)&time_vec_full==time_index(t)) = weight_array(t,a);
        end
    end
    wt_vec_full(1~=flag_vec_new) = NaN;
    regression_set.([dp_var_flag '_wt']) = wt_vec_full;          
end
% Save
writetable(regression_set, [inputOutputPath 'final_regression_set.csv'])