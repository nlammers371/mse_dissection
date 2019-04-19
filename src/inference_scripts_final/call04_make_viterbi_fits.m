% Script to generate Viterbi Fits for Inference traces
close all
clear 
addpath('..\utilities\');

project = 'mHMMeve2_weka_inf_2018_05_07'; %project identifier
figure_id = 'main_figure5';
w = 7; %memory assumed for inference
K = 3; %states used for final inference
Tres = 20; %time resolution

infPath = ['../../out/' figure_id '/' project '/'];
tracePath = ['../../dat/figure_data/' project '/'];
%%% Load Data
load([infPath '\hmm_results_summary.mat'])
load([tracePath '/inference_traces_' project '_dT' num2str(Tres) '.mat'])


%%% Make Viterbi Fits
hmm_regions = [hmm_results.binID];
alpha = hmm_results(1).alpha;
dT = hmm_results(1).dT;
viterbi_fit_struct = struct;
skipped_stripes = [];

parfor i = 1:length(trace_struct_final)
    MeanAP = round(trace_struct_final(i).MeanAP);        
    % if trace is from far anterior/posterior flank, use results for
    % closest AP bin
    [~,hmm_index] = min(abs(hmm_regions-MeanAP));    
    hmm_bin = hmm_results(hmm_index);

    viterbi_fit_struct(i).skipped = 0;
    viterbi_fit_struct(i).hmm_bin = hmm_regions(hmm_index);
    if isempty(hmm_bin)|| length(trace_struct_final(i).fluo_interp) < w        
        viterbi_fit_struct(i).skipped = 1;
        viterbi_fit_struct(i).ParticleID = trace_struct_final(i).ParticleID;
        continue
    end    
    v = hmm_bin.initiation_mean/60*dT;
    noise = hmm_bin.noise_mean;
    pi0_log = log([1;1;1]/3);
    A_log = reshape(log(hmm_bin.A_mean),K,K);                
    fluo = trace_struct_final(i).fluo_interp;            
    v_fit = viterbi (fluo, v', noise, pi0_log,A_log, K, w, alpha);    
    v_fit.time_exp = trace_struct_final(i).time_interp;
    v_fit.fluo_exp = fluo;            
    v_fit.v = v;
    v_fit.w = w;        
    v_fit.alpha = alpha;
    v_fit.noise = noise;
    v_fit.pi0 = exp(pi0_log);
    v_fit.A = exp(A_log);
    viterbi_fit_struct(i).v_fit = v_fit;
    viterbi_fit_struct(i).ParticleID = trace_struct_final(i).ParticleID;
    disp([num2str(i) ' of ' num2str(length(trace_struct_final)) ' completed'])
end

save([infPath '/viterbi_fits.mat'],'viterbi_fit_struct') 