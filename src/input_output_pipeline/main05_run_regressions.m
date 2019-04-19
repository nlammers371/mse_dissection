% Script to systematically fit nonlinear regression models to experimental
% data
function [use_weights,constrained_inference,ap_raw_flag] = main05_run_regressions(project,inference_id,varargin)

addpath('../utilities');
% data paths
inputOutputPath = ['../../out/input_output/' project '/'];
%%%% Set Fit Parameters
% indicate whether to use raw or adjusted AP coordinates 
ap_raw_flag = 0;
% define response variables to test
dp_variables = {'on','z_on'};
dp_var_flags = {'bin_act_flag','z_on_flag'};
% Specify list of predictors to use
pd_variables = {'Bcd','Hb','Gt','Kr'};
% define max power of terms allowed (eg 2=[Hb]^2)
max_power_vec = [1 1 1 1]; %[6 6 6 6]; 
% max number of times variations of same TF may appear in regression
max_multi_vec = [1 1 1 1]; 
% bootstrap params
bootstrap_flag = 1;
n_boots = 100;
boot_samp_size = 50000;
% cross-valiudation params
cv_flag = 1;
n_cv_groups = 10;
% specify whether or not to use sampling weights
use_weights = 1;
% set upper and lower bounds
max_vars = 1+numel(pd_variables);
min_vars = 1;


for i = 1:numel(varargin)
    if ischar(varargin{i})
        eval([varargin{i} ' = varargin{i+1};'])
    end
end


% load data
regression_set = readtable([inputOutputPath 'final_regression_set.csv']);
load([inputOutputPath 'tf_input_struct.mat']);
% specify constraints on tf slope params
if strcmpi(inference_id, 'main_text')
    constrained_inference = 1;
    ub_base = [Inf Inf 0 0];
    lb_base = [0 0 -Inf -Inf];
elseif strcmpi(inference_id, 'appendix')    
    constrained_inference = 0;
    ub_base = [Inf Inf Inf Inf];
    lb_base = [-Inf -Inf -Inf -Inf];
else 
    warning('Inference ID not recognized. No constraints applied')
    constrained_inference = 0;
    ub_base = [Inf Inf Inf Inf];
    lb_base = [-Inf -Inf -Inf -Inf];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add TF data to eve table
ap_field_eve2 = 'AP';
if ap_raw_flag
    ap_field_eve2 = 'APOrig';
end
InterpGrid = tf_input_struct(1).InterpGrid;
ap_vec = tf_input_struct(1).ap_vec;
final_regression_set = add_tf_fields(regression_set,tf_input_struct,ap_field_eve2);
% generate additional fields for higher-order terms as needed
reg_variables = {};%pd_variables;
ub_full = [];
lb_full = [];
for i = 1:numel(pd_variables)
    tf_vec = final_regression_set.(pd_variables{i});
    for j = 1:max_power_vec(i)      
        tf_str = [pd_variables{i} '_' num2str(j)];
        reg_variables = [reg_variables{:} {tf_str}];                      
        ub_full = [ub_full ub_base(i)];
        lb_full = [lb_full lb_base(i)];
    end
end
final_regression_set = augment_structure(final_regression_set,pd_variables,reg_variables);

% generate project string
inference_string = [inference_id '_wt' num2str(use_weights) '_apRaw' num2str(ap_raw_flag) '_ct' num2str(constrained_inference)];
regressionPath = [inputOutputPath '/results_' inference_string '/'];
mkdir(regressionPath)
% ref indice
reg_indices = 1:numel(reg_variables);

% queue up all variable combinations that adhere to specified constraints
reg_var_cell = {};
reg_ind_cell = {};
for n = min_vars:max_vars
    combos = combnk(reg_indices,n);        
    for i = 1:size(combos,1)
        n_vec = [];
        for k = 1:numel(pd_variables)
            n_vec = [n_vec sum(contains(reg_variables(combos(i,:)),pd_variables{k}))];
        end
        if sum(n_vec>max_multi_vec)==0
            reg_var_cell = [reg_var_cell {reg_variables(combos(i,:))}];
            reg_ind_cell = [reg_ind_cell{:} {combos(i,:)}];
        end
    end    
end

% now iterate through response variables and model architectures
reg_struct = struct;
for rg = 1:numel(dp_variables) % variables
    tic
    % extract variable info
    dp_var = dp_variables{rg};
    dp_var_flag = dp_var_flags{rg};    
    reg_struct(rg).dp_var = dp_var;
    reg_struct(rg).dp_flag = dp_var_flag;   
    % extract weights
    dp_wt_vec = final_regression_set.([dp_var_flag '_wt']);    
    reg_struct(rg).base_variables = pd_variables; 
    reg_struct(rg).reg_variables = reg_variables;           
    reg_struct(rg).max_power_vec = max_power_vec;
    reg_struct(rg).max_multi_vec = max_multi_vec;
    reg_struct(rg).ap_field_eve2 = ap_field_eve2;
    reg_struct(rg).use_weights = use_weights;
    reg_struct(rg).sample_size = boot_samp_size;
    temp_struct = struct;
    % iterate through models
    for md = 1:numel(reg_var_cell) 
        pd_vars = reg_var_cell{md};
        reg_indices = reg_ind_cell{md};
        Y_flag = final_regression_set.(dp_var_flag);        
        % response var
        Y = final_regression_set.(dp_var); 
        set_vec = final_regression_set.Set;
        % get predictors
        X = []; 
        for n = 1:numel(pd_vars)
            X = [X final_regression_set.(pd_vars{n})];
        end
        % remove NaNs
        X_flag = sum(isnan(X),2)==0;
        Y = Y(Y_flag==1&X_flag==1);
        X = X(Y_flag==1&X_flag==1,:);  
        set_vec = set_vec(Y_flag==1&X_flag==1); 
        set_index = unique(set_vec);
        var_wt_vec = dp_wt_vec(Y_flag==1&X_flag==1,:);
        var_wt_vec(isnan(var_wt_vec)) = 0;
        
        n_sample_points = min(boot_samp_size,numel(Y));        
        set_resamp = true;
        if ~bootstrap_flag
%             n_sample_points = numel(Y);
            n_boots = 1;
            set_resamp = false;
        end
        % hold sample size constant across bootstraps
        % sampling index vector    
        beta_mat = NaN(size(X,2)+1,n_boots);
        dev_vec_train = NaN(1,n_boots);     
        dev_vec_test = NaN(1,n_boots);     

        for n = 1:n_boots
            % first bootstrap at level of sets
            % training data
            set_boot_vec = randsample(set_index,numel(set_index),set_resamp);
            boot_index_vec = []; 
            for s = 1:numel(set_boot_vec)
                boot_index_vec = [boot_index_vec find(set_vec==set_boot_vec(s))'];
            end
            boot_wt_vec = var_wt_vec(boot_index_vec);
            % draw sample
            if use_weights == 1
                sample_ids = randsample(boot_index_vec,n_sample_points,true,boot_wt_vec);
            else
                sample_ids = randsample(boot_index_vec,n_sample_points,true);
            end
            X_boot = X(sample_ids,:);
            Y_boot = Y(sample_ids,:);
            if cv_flag                
                cv_samp_size = ceil(n_sample_points/n_cv_groups);                
                % randomly assign data rows to groups
                cv_index = 1:n_cv_groups;
                cv_vec_raw = randsample(repelem(cv_index,cv_samp_size),numel(repelem(cv_index,cv_samp_size)),false);
                cv_vec = cv_vec_raw(1:n_sample_points);
            else
                cv_vec = zeros(1,n_sample_points);
                n_cv_groups = 1;
            end
            cv_beta_mat = NaN(size(X,2)+1,n_cv_groups);
            cv_dev_vec_train = NaN(1,n_cv_groups);     
            cv_dev_vec_test = NaN(1,n_cv_groups);  
            for cv = 1:n_cv_groups
                train_ft = cv_vec~=cv;
                test_ft = cv_vec==cv;
                 
                X_train = X_boot(train_ft,:);
                Y_train = Y_boot(train_ft,:); 
                X_test = X_boot(test_ft,:);
                Y_test = Y_boot(test_ft,:); 
                X_train_ct = [ones(size(X_train,1),1) X_train];                
                 
                fcn = @(b) objfunc(b',X_train_ct,Y_train);                                   

                lb = [-Inf lb_full(reg_indices)];
                ub = [Inf ub_full(reg_indices)];   

                x0 = rand(1,numel(reg_indices)+1).*(sign(lb)+sign(ub))...
                    -rand(1,numel(reg_indices)+1).*~(sign(lb)+sign(ub))*10;   
                options = optimoptions('fmincon','Display','off');
                [B, ~] = fmincon(fcn,x0,[],[],[],[],lb,ub,[],options);%X_train_ct,double(Y_train));                                                                                   
                B = -B';      
                % evaluate model in context of 
                pi_test = mnrval(B,X_test);    
                cv_dev_vec_test(cv) = mean([log(pi_test(Y_test==0,1))' log(pi_test(Y_test==1,2))']);  
                % evaluate model performance for training data
                pi_train = mnrval(B,X_train);    
                cv_dev_vec_train(cv) = mean([log(pi_train(Y_train==0,1))' log(pi_train(Y_train==1,2))']);  
                % record coefficients
                cv_beta_mat(:,cv) = B;                                                                                         
            end   
            beta_mat(:,n) = nanmean(cv_beta_mat,2);
            dev_vec_train(n) = nanmean(cv_dev_vec_train);
            dev_vec_test(n) = nanmean(cv_dev_vec_test);
        end
        % record info
        temp_struct(md).pd_vars = pd_vars;        
        temp_struct(md).n_vars = numel(pd_vars);        
        temp_struct(md).pd_indices = reg_indices;
        temp_struct(md).beta_fit = nanmean(beta_mat,2) ;                
        temp_struct(md).beta_mat = beta_mat;
        temp_struct(md).n_samples = n_sample_points;
        temp_struct(md).score_train = nanmean(dev_vec_train)/2; 
        temp_struct(md).score_test = nanmean(dev_vec_test)/2; 
        temp_struct(md).score_vec_train = dev_vec_train/2; 
        temp_struct(md).score_vec_test = dev_vec_test/2; 
        temp_struct(md).constrained = constrained_inference;
        % sampling info
        temp_struct(md).cv_flag = cv_flag;
        temp_struct(md).bootstrap_flag = bootstrap_flag;
        % constraint info      
        temp_struct(md).ub = ub_full(reg_indices);
        temp_struct(md).lb = lb_full(reg_indices);                            
    end
    reg_struct(rg).reg_results = temp_struct;
    disp(['Completed Regressions for ' dp_var_flag ' (' num2str(round(toc)) ' sec)'])    
end
% save

save([regressionPath 'input_output_results.mat'],'reg_struct')