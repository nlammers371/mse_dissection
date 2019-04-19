% Compile results into summary structure. Generate summary plots
function compile_inference_results(project,inference_type,figure_id,K,w)
addpath('../utilities/');
%------------------------------Set System Params--------------------------%
% set defaults
Tres = 20; %Time Resolution
mindp = 500; % minimum # data points permitted for use
ap_groups = {26:29,30:33,34:37,38:39,40:41,42:43,44:47,48:52};
%-----------------------------ID Variables--------------------------------%
%Generate filenames and writepath
file_id_string =  ['eveSet_w' num2str(w) '*K' num2str(K) '*']; 
dataPath = ['../../out/' inference_type '/' project '/' figure_id '/'];
outPath = ['../../out/' inference_type '/' project '/'];
mkdir(outPath)

%---------------------------------Read in Files---------------------------%
files = dir([dataPath file_id_string]);
filenames = {};
for i = 1:length(files)
    filenames = [filenames {files(i).name}]; 
end

if isempty(filenames)
    warning('No file with specified inference parameters found')
    return
end

%Iterate through result sets and combine into 1 master struct
inference_results = struct;
f_pass = 1;

for f = 1:length(filenames) 
    load([dataPath filenames{f}]);

    if length(fieldnames(output)) < 2
        continue
    elseif output.N < mindp
        continue
%         newDir = [dataPath './insufficient_data/'];
%         mkdir(newDir)
%         fp = [dataPath filenames{f}];
%         movefile(fp, newDir)
    end
    ap_flag = false;
    for i = 1:numel(ap_groups)
        ap_flag = ap_flag || all(ismember(output.APbin,ap_groups{i}));
    end
    if ap_flag
        for fn = fieldnames(output)'
            inference_results(f_pass).(fn{1}) = output.(fn{1});
        end
        inference_results(f_pass).source = filenames{f};        
        f_pass = f_pass + 1
    end
end

% Calculate summary statistics
% Adjust rates as needed (address negative off-diagonal values)
% Define convenience Arrays and vectors
alpha = inference_results(1).alpha;
% remove entries with too few data points
n_dp_vec = [inference_results.N];
inference_results = inference_results(n_dp_vec>=mindp);

% get position group info
bin_vec = [];
bin_cell = {};
for i = 1:length(inference_results)
    bin_vec = [bin_vec mean(inference_results(i).APbin)];
    bin_cell = [bin_cell{:} {inference_results(i).APbin}];
end
bin_range = unique(bin_vec);

% time info
time_vec = [inference_results.t_inf];
time_index = unique(time_vec);

% initialize arrays to store inference parameters
initiation_rates = NaN(K,length(inference_results)); % r
pi0_all = NaN(K,length(inference_results));    
noise_all = NaN(1,length(inference_results)); % sigma  
n_dp_all = NaN(1,length(inference_results));
n_traces_all = NaN(1,length(inference_results));
effective_times_all = NaN(1,length(inference_results));
dwell_all = NaN(K,length(inference_results));

% Extract variables and perform rate fitting as needed
for i = 1:length(inference_results)    
    [initiation_rates(:,i), ranked_r] = sort(60*[inference_results(i).r]); 
    pi0 = inference_results(i).pi0(ranked_r);    
    noise_all(i) = sqrt(inference_results(i).noise); 
    n_dp_all(i) = inference_results(i).N;
    n_traces_all(i) = length(inference_results(i).particle_ids);
    effective_times_all(i) = nanmean([inference_results(i).time_data{:}]);
    A = reshape(inference_results(i).A,K,K);
    A = A(ranked_r, ranked_r);
    %Obtain raw R matrix
    R = prob_to_rate(A,Tres);
    inference_results(i).R_mat = R;
    Rcol = reshape(R,1,[]);
    inference_results(i).R = Rcol;
    %Check for imaginary and negative elements. If present, perform rate
    %fit   
    inference_results(i).r_fit_flag = 0;    
    if ~isreal(Rcol) || sum(Rcol<0)>K
        inference_results(i).r_fit_flag = 1;
        out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 1);            
        inference_results(i).R_fit = 60*out.R_out;
        r_diag = 60*diag(out.R_out);
    elseif false %sum(Rcol<0)>K
        glb_all(i).r_fit_flag = 2;
        R_conv = 60*(R - eye(K).*R);
        R_conv(R_conv<0) = 0;
        R_conv(eye(K)==1) = -sum(R_conv);
        glb_all(i).R_fit = R_conv;
        r_diag = diag(R_conv);
    else
        inference_results(i).R_fit = 60*inference_results(i).R_mat;
        r_diag = 60*diag(inference_results(i).R_mat); 
    end    
    dwell_all(:,i) = -1./r_diag;
end

% Save rate info
%Make 3D arrays to stack matrices
R_orig_array = NaN(length(inference_results),K^2);
R_fit_array = NaN(length(inference_results),K^2);
A_array = NaN(length(inference_results),K^2);
for i = 1:length(inference_results)
    R_fit_array(i,:) = reshape(inference_results(i).R_fit,1,[]);
    R_orig_array(i,:) = reshape(real(inference_results(i).R_mat),1,[]);
    A_array(i,:) = reshape(real(inference_results(i).A_mat),1,[]);
end

% 2D Arrays to store moments (A's are calculated to have for output structure)
med_R_orig = zeros(length(time_index),K^2,length(bin_range));
std_R_orig = zeros(length(time_index),K^2,length(bin_range));
avg_A = zeros(length(time_index),K^2,length(bin_range));
std_A = zeros(length(time_index),K^2,length(bin_range));
med_R_fit = zeros(length(time_index),K^2,length(bin_range));
std_R_fit = zeros(length(time_index),K^2,length(bin_range));
med_on_ratio = NaN(length(time_index),length(bin_range));
std_on_ratio = NaN(length(time_index),length(bin_range));
med_off_ratio = NaN(length(time_index),length(bin_range));
std_off_ratio = NaN(length(time_index),length(bin_range));
for b = 1:length(bin_range)
    bin = bin_range(b);
    for t = 1:length(time_index)
        time = time_index(t);
        A_mean = mean(A_array(bin_vec==bin&time_vec==time,:),1);
        A_reshape = reshape(A_mean,K,K);
        
        % normalize A
        A_mean = reshape(A_reshape ./ repmat(sum(A_reshape),K,1),1,[]);
        avg_A(t,:,b) = A_mean;   
        
        %Calculate and Store R moments        
                med_R_orig(t,:,b) = median(R_orig_array(bin_vec==bin&time_vec==time,:),1);        
        std_R_orig(t,:,b) = .5*(quantile(R_orig_array(bin_vec==bin&time_vec==time,:),.75,1)-...
                                quantile(R_orig_array(bin_vec==bin&time_vec==time,:),.25,1));            
        med_R_fit(t,:,b) = median(R_fit_array(bin_vec==bin&time_vec==time,:),1);
        std_R_fit(t,:,b) = .5*(quantile(R_fit_array(bin_vec==bin&time_vec==time,:),.75,1)-...
                               quantile(R_fit_array(bin_vec==bin&time_vec==time,:),.25,1));
        if K == 3
            k_on_rat = R_fit_array(bin_vec==bin&time_vec==time,2)./...
                            R_fit_array(bin_vec==bin&time_vec==time,6);
            med_on_ratio(t,b) = median(k_on_rat);
            std_on_ratio(t,b) = .5*(quantile(k_on_rat,.75,1)-...
                                   quantile(k_on_rat,.25,1));
            k_off_rat = R_fit_array(bin_vec==bin&time_vec==time,8)./...
                            R_fit_array(bin_vec==bin&time_vec==time,4);
            med_off_ratio(t,b) = median(k_off_rat);
            std_off_ratio(t,b) = .5*(quantile(k_off_rat,.75,1)-...
                                   quantile(k_off_rat,.25,1));
        end
    end
end

% Calculate state occupancy 
occupancy = zeros(K,length(inference_results));
for i = 1:length(inference_results)
    [~, ranked_r] = sort([inference_results(i).r]);
    A = inference_results(i).A_mat;
    A = A(ranked_r,ranked_r);
    [V,D] = eig(A);
    ind = diag(D)==max(diag(D));
    steady = V(:,ind)./sum(V(:,ind));
    occupancy(:,i) = steady;
    inference_results(i).occupancy = steady;
end

% initialize additional parameter arrays
med_dwell = NaN(K,length(bin_range),length(time_index));
std_dwell = NaN(K,length(bin_range),length(time_index));
med_initiation = NaN(K,length(bin_range),length(time_index));
std_initiation = NaN(K,length(bin_range),length(time_index));
med_init_ratio = NaN(length(time_index),length(bin_range));
std_init_ratio = NaN(length(time_index),length(bin_range));
med_occupancy = NaN(K,length(bin_range),length(time_index));
std_occupancy = NaN(K,length(bin_range),length(time_index));
avg_pi0 = NaN(K,length(bin_range),length(time_index));
std_pi0 = NaN(K,length(bin_range),length(time_index));
med_noise = NaN(1,length(bin_range),length(time_index));
std_noise = NaN(1,length(bin_range),length(time_index));
med_n_dp = NaN(1,length(bin_range),length(time_index));
med_n_tr = NaN(1,length(bin_range),length(time_index));
med_eff_time = NaN(1,length(bin_range),length(time_index));
n_boots_total = NaN(length(time_index),length(bin_range));
for i = 1:length(bin_range)
    bin = bin_range(i);
    for t = 1:length(time_index)
        time = time_index(t);
        n_boots_total(t,i) = sum(bin_vec==bin&time_vec==time);
        for k = 1:K
            if ~isempty(initiation_rates(k,bin_vec==bin&time_vec==time))
                med_initiation(k,i,t) =  median(initiation_rates(k,bin_vec==bin&time_vec==time));            
                std_initiation(k,i,t) = .5*(quantile(initiation_rates(k,bin_vec==bin&time_vec==time),.75)-...
                        quantile(initiation_rates(k,bin_vec==bin&time_vec==time),.25));                  
                med_occupancy(k,i,t) = median(occupancy(k,bin_vec==bin&time_vec==time));
                std_occupancy(k,i,t) = .5*(quantile(occupancy(k,bin_vec==bin&time_vec==time),.75)-...
                        quantile(occupancy(k,bin_vec==bin&time_vec==time),.25)); 
                med_dwell(k,i,t) = median(dwell_all(k,bin_vec==bin&time_vec==time));
                std_dwell(k,i,t) = .5*(quantile(dwell_all(k,bin_vec==bin&time_vec==time),.75)-...
                        quantile(dwell_all(k,bin_vec==bin&time_vec==time),.25));        
            end
        end     
        if K == 3
            init_rat = initiation_rates(3,bin_vec==bin&time_vec==time)./...
                            initiation_rates(2,bin_vec==bin&time_vec==time);
            med_init_ratio(t,i) = median(init_rat);
            std_init_ratio(t,i) = .5*(quantile(init_rat,.75)-...
                        quantile(init_rat,.25));
        end
        med_noise(i,t) = median(noise_all(bin_vec==bin&time_vec==time));
        std_noise(i,t) = .5*(quantile(noise_all(bin_vec==bin&time_vec==time),.75)-...
                    quantile(noise_all(bin_vec==bin&time_vec==time),.25));
        med_n_dp(i,t) = median(n_dp_all(bin_vec==bin&time_vec==time));
        med_n_tr(i,t) = median(n_traces_all(bin_vec==bin&time_vec==time));
        med_eff_time(i,t) = median(effective_times_all(bin_vec==bin&time_vec==time));
    end
end

% effective off rate, on rate and emission rates 
if K == 3
    on_rate_eff_med = NaN(length(time_index),length(bin_range));
    on_rate_eff_ste = NaN(length(time_index),length(bin_range));
    off_rate_eff_med = NaN(length(time_index),length(bin_range));
    off_rate_eff_ste = NaN(length(time_index),length(bin_range));
    init_rate_eff_med = NaN(length(time_index),length(bin_range));
    init_rate_eff_ste = NaN(length(time_index),length(bin_range));
    on_occ_eff_med = NaN(length(time_index),length(bin_range));
    on_occ_eff_ste = NaN(length(time_index),length(bin_range));
    for b = 1:length(bin_range)
        bin = bin_range(b);
        for t = 1:length(time_index)
            time = time_index(t);        
            %Calculate and Store R moments  
            r21_vec = R_fit_array(bin_vec==bin&time_vec==time,2);
            r31_vec = R_fit_array(bin_vec==bin&time_vec==time,3);
            r12_vec = R_fit_array(bin_vec==bin&time_vec==time,4);
%             r32_vec = R_fit_array(bin_vec==bin&time_vec==time,6);
            r13_vec = R_fit_array(bin_vec==bin&time_vec==time,7);            
%             r23_vec = R_fit_array(bin_vec==bin&time_vec==time,8)/2;
            pi1_vec = occupancy(1,bin_vec==bin&time_vec==time)';
            pi2_vec = occupancy(2,bin_vec==bin&time_vec==time)';
            pi3_vec = occupancy(3,bin_vec==bin&time_vec==time)';
            % occupancy first
            pi_on_all = pi2_vec + pi3_vec;
            on_occ_eff_med(t,b) = median(pi_on_all);
            on_occ_eff_ste(t,b) = .5*(quantile(pi_on_all,.75)-quantile(pi_on_all,.25));
            % transition rates
            on_rate_all = r21_vec + r31_vec;
            on_rate_eff_med(t,b) = median(on_rate_all);
            on_rate_eff_ste(t,b) = .5*(quantile(on_rate_all,.75)-quantile(on_rate_all,.25));
            off_rate_all = (pi2_vec.*r12_vec+pi3_vec.*r13_vec) ./ (pi2_vec+pi3_vec);
            off_rate_eff_med(t,b) = median(off_rate_all);
            off_rate_eff_ste(t,b) = .5*(quantile(off_rate_all,.75)-quantile(off_rate_all,.25));
            init_rate_all =  (pi2_vec.*initiation_rates(2,bin_vec==bin&time_vec==time)'+...
                pi3_vec.*initiation_rates(3,bin_vec==bin&time_vec==time)')./(pi2_vec+pi3_vec);  
            init_rate_eff_med(t,b) = median(init_rate_all);
            init_rate_eff_ste(t,b) = .5*(quantile(init_rate_all,.75)-quantile(init_rate_all,.25));
        end
    end
end
%%% Make Output Struct With Relevant Fields
hmm_results = struct;
for i = 1:length(bin_range)
    for t = 1:length(time_index)
        ind = (i-1)*length(time_index) + t;
        hmm_results(ind).n_boots_total = n_boots_total(t,i);
        hmm_results(ind).initiation_mean = med_initiation(:,i,t);
        hmm_results(ind).initiation_std = std_initiation(:,i,t);    
        hmm_results(ind).occupancy_mean = med_occupancy(:,i,t);
        hmm_results(ind).occupancy_std = std_occupancy(:,i,t);        
        hmm_results(ind).pi0_mean = avg_pi0(:,i,t);
        hmm_results(ind).pi0_std = std_pi0(:,i,t);
        hmm_results(ind).dwell_mean = med_dwell(:,i,t);
        hmm_results(ind).dwell_std = std_dwell(:,i,t);    
        hmm_results(ind).A_mean = avg_A(t,:,i);        
        hmm_results(ind).R_orig_mean = med_R_orig(t,:,i);
        hmm_results(ind).R_orig_std = std_R_orig(t,:,i);
        hmm_results(ind).R_fit_mean = med_R_fit(t,:,i);
        hmm_results(ind).R_fit_std = std_R_fit(t,:,i);
        hmm_results(ind).noise_mean = med_noise(i,t);
        hmm_results(ind).noise_std = std_noise(i,t);
        hmm_results(ind).N_dp = med_n_dp(i,t);
        hmm_results(ind).N_tr = med_n_tr(i,t);
        hmm_results(ind).t_inf_effective = med_eff_time(i,t);
        hmm_results(ind).noise_mean = med_noise(i,t);
        if K == 3
            hmm_results(ind).on_occ_eff_mean = on_occ_eff_med(t,i);
            hmm_results(ind).on_occ_eff_ste = on_occ_eff_ste(t,i);
            hmm_results(ind).on_ratio_mean = med_on_ratio(t,i);
            hmm_results(ind).on_ratio_ste = std_on_ratio(t,i);
            hmm_results(ind).off_ratio_mean = med_off_ratio(t,i);
            hmm_results(ind).off_ratio_ste = std_off_ratio(t,i);
            hmm_results(ind).init_ratio_mean = med_init_ratio(t,i);
            hmm_results(ind).init_ratio_ste = std_init_ratio(t,i);
            hmm_results(ind).eff_on_rate = on_rate_eff_med(t,i);
            hmm_results(ind).eff_on_ste = on_rate_eff_ste(t,i);
            hmm_results(ind).eff_off_rate = off_rate_eff_med(t,i);
            hmm_results(ind).eff_off_ste = off_rate_eff_ste(t,i);
            hmm_results(ind).eff_init_rate = init_rate_eff_med(t,i);
            hmm_results(ind).eff_init_ste = init_rate_eff_ste(t,i);
        end
        hmm_results(ind).binID = bin_range(i);
        hmm_results(ind).t_inf = time_index(t);
        hmm_results(ind).alpha = alpha;
        hmm_results(ind).dT = Tres;
        hmm_results(ind).w = w;
    end
end
save([outPath 'hmm_results_summary_' figure_id '.mat'],'hmm_results')
if K == 3 && w == 7
    outPath2 = ['../../out/appendix8_figure2w7/' project '/'];
    mkdir(outPath2);
    save([outPath2 'hmm_results_summary_' figure_id '.mat'],'hmm_results')
end
% NL: these files get too large, not saving
% save([outPath 'hmm_results_compiled.mat'],'inference_results')
                 