% Script to assess relative impact of parameters for each regression model
% and flag general "motifs" of interest
function main06_add_trace_model_motifs(project,inference_id,use_weights,constrained_inference,ap_raw_flag)

% data paths
inputOutputPath = ['../../out/input_output/' project '/'];
inference_string = [inference_id '_wt' num2str(use_weights) '_apRaw' num2str(ap_raw_flag) '_ct' num2str(constrained_inference)];
regressionPath = [inputOutputPath '/results_' inference_string '/'];
% load regression data
regression_set = readtable([inputOutputPath 'final_regression_set.csv']);
load([inputOutputPath 'tf_input_struct.mat']);

% load reg results structure
load([regressionPath 'input_output_results.mat']);
reg_struct_out = reg_struct;
% build full regression set
ap_field_eve2 = reg_struct(1).ap_field_eve2;


temp_regression_set = add_tf_fields(regression_set,tf_input_struct,ap_field_eve2);
% add fields on the fly as needed
reg_variables = reg_struct(1).reg_variables;
base_variables = reg_struct(1).base_variables;
final_regression_set = augment_structure(temp_regression_set,base_variables,reg_variables);

% extract useful variables
n_var_vec = [reg_struct(1).reg_results.n_vars];

% define a set of "motifs" that are of interest
motif_string_cell = {{'Gt','Kr'},{'Gt','Kr','Hb'},{'Gt','Kr','Bcd'},{'Gt','Kr','Hb'},...
    {'Gt','Kr','Bcd','Bcd'},{'Gt','Kr','Hb','Hb'},{'Gt','Kr','Bcd','Bcd','Hb'}};
motif_sign_cell = {[-1 -1], [-1 -1 1], [-1 -1 -1], [-1 -1 -1], [-1 -1 0 0], [-1 -1 0 0],...
    [-1 -1 0 0 1]};
motif_name_cell = {'dual repressor', 'dual repressor Hb activate', 'dual repressor Bcd repress',...
    'dual repressor Hb repress', 'dual repressor bifunctional Bcd', 'dual repressor bifunctional Hb','dual repressor bifunctional Bcd Hb act'};
motif_index_vec = 1:numel(motif_string_cell);

% iterate through all models to assess contribution from each predictor
for i = 1:numel(reg_struct)
    dp_var = reg_struct(i).dp_var;    
    dp_flag = reg_struct(i).dp_flag;        
    dp_wt_vec = final_regression_set.([dp_flag '_wt']);
    reg_results = reg_struct(i).reg_results;    
    use_weights = reg_struct(i).use_weights;       
    for j = 1:numel(reg_results)
        pd_indices = reg_results(j).pd_indices;
        pd_vars = reg_results(j).pd_vars;
        % extract regression arrays
        Y = final_regression_set.(dp_var);  
        
        AP = final_regression_set.AP;  
        Time = final_regression_set.Time;  
        
        Y_flag = final_regression_set.(dp_flag);  
        % get predictors
        X = []; 
        for n = 1:numel(pd_vars)
            X = [X final_regression_set.(pd_vars{n})];
        end
        % remove NaNs
        X_flag = sum(isnan(X),2)==0;
        Y = Y(Y_flag==1&X_flag==1);
        X = X(Y_flag==1&X_flag==1,:);        
        var_wt_vec = dp_wt_vec(Y_flag==1&X_flag==1,:);
        var_wt_vec(isnan(var_wt_vec)) = 0;
        AP = AP(Y_flag==1&X_flag==1,:);
        Time = Time(Y_flag==1&X_flag==1,:);
        % extract reg result info
        
        beta = reg_results(j).beta_fit;        
        sign_vec = sign(-beta');        
        cb_vec = NaN(size(pd_indices));
        b_index = 0:numel(pd_indices);
        % assess penalty of removing each variable in turn
        ns = min(reg_struct(i).sample_size,size(X,1));
        index_vec= 1:numel(Y);
        if use_weights == 1
            test_vec = randsample(index_vec,ns,true,var_wt_vec);
        else
            test_vec = randsample(index_vec,ns,false);
        end
        X_test = X(test_vec,:);
        Y_test = Y(test_vec,:);
        pi_test = mnrval(beta,X_test);    
        pi_vec = [log(pi_test(Y_test==0,1))' log(pi_test(Y_test==1,2))'];
        base_score = sum(pi_vec);         
        % let's get error-weighted position and times while we're at it
        ap_test = AP(test_vec);
        time_test = Time(test_vec);        
        ap_err = sum((pi_vec.*[ap_test(Y_test==0)' ap_test(Y_test==1)'])) / base_score;
        time_err = sum((pi_vec.*[time_test(Y_test==0)' time_test(Y_test==1)'])) / base_score;         
        
        for k = 1:numel(pd_indices)            
            % obtain predictions
            X_sub = X_test(:,b_index(2:end)~=k);
            beta_sub = beta(b_index~=k);
            % catch single predictor cases
            if isempty(X_sub)
                X_sub = zeros(size(Y_test));
                beta_sub = [beta(1) 0]';
            end
            pi_sub = mnrval(beta_sub,X_sub);    
            sub_dev = sum([log(pi_sub(Y_test==0,1))' log(pi_sub(Y_test==1,2))']);
            cb_vec(k) = round(base_score - sub_dev,1);            
        end        
        % convert to percent...should I exponentiate?
        cb_vec_norm = cb_vec;
        cb_vec_norm(cb_vec_norm<0) = 0;
        cb_vec_norm = cb_vec_norm / sum(cb_vec_norm);
        cb_vec_norm(cb_vec==0) = 0;
        reg_results(j).pd_contributions_raw = [NaN cb_vec];
        reg_results(j).pd_contributions_norm = [NaN cb_vec_norm];
        reg_results(j).sign_vec = sign_vec;
        % assign motif (if relevant)
        % motif match declared if set of motif variables accounts for >95%
        % of model score and there is not a smaller subset that surpasses
        % that threshold        
        motif_index = 0;        
        sign_vec = sign_vec(2:end);        
        for k = motif_index_vec
            pd_var_sub = pd_vars;
            motif_strings = motif_string_cell{k};
            motif_signs = motif_sign_cell{k};
            sig_vec = [];
            for l = 1:numel(motif_strings)
                ct = contains(pd_var_sub,motif_strings{l});                
                if sum(ct)>0&&motif_signs(l)==sum(sign_vec(ct))
                    sig_vec = [sig_vec cb_vec_norm(ct)];                   
                    pd_var_sub(ct) = {'NA'};
                end
            end            
            % check to
            if numel(sig_vec) == numel(motif_strings)
                sig = sum(sig_vec);
                sig_less = sig - min(sig_vec);                
                if (sig_less <.95 && sig >= .95) && motif_index == 0
                    motif_index = k;
                elseif (sig_less <.95 && sig >= .95) && motif_index ~= 0
                    error('degenerate motif assignment')
                end
            end
        end            
        reg_results(j).ap_err = ap_err;        
        reg_results(j).time_err = time_err;        
        reg_results(j).motif_index = motif_index;        
        reg_results(j).sig_vec = cb_vec_norm>.05;        
    end    
    reg_struct_out(i).reg_results = reg_results;
    reg_struct_out(i).motif_string_cell = motif_string_cell;
    reg_struct_out(i).motif_sign_cell = motif_sign_cell;
    reg_struct_out(i).motif_name_cell = motif_name_cell;            
end
reg_struct = reg_struct_out;
save([regressionPath '/input_output_results_final.mat'], 'reg_struct')