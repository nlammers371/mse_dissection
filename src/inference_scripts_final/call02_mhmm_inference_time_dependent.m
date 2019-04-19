% Master script to call mHMM inference for various conditions of interest
clear
close all
% First define variables that are common across all inference runs
project = 'mHMMeve2_weka_inf_2018_05_07'; % Project identifier
Tres = 20; % Time Resolution
t_window = 15;
inference_times = (7.5:2.5:45);
% inference_times = randsample(inference_times,numel(inference_times),false);
savio = 0;
if savio
    outPath = '/global/scratch/nlammers/temporal_inference_results/'; 
else    
    outPath = '../../out/temporal_inference_results/';
end
mkdir(outPath);
% load inference traces 
dataPath = ['../../dat/figure_data/' project '/']; %Path to raw data
dataName = ['inference_traces_' project '_dT' num2str(Tres) '.mat'];
load([dataPath dataName]);
% designate inference regions
ap_ref_cell = {-1:1};%{-7:-4,-3:-2,-1:1,2:3,4:7}; 
bin_groups = ap_ref_cell;
if savio
    %Get environment variable from job script
    savio_groups = {str2num(getenv('SLURM_ARRAY_TASK_ID'))};    
    bin_groups = cell(1,length(savio_groups));
    for i =1:length(bin_groups)
        bin_groups{i} = ap_ref_cell{savio_groups{i}};
    end
end  
 
% %Figure 6: 3 State time-varying inference 
% w = 7;
% K = 3;
% id_string = 'main_figure6';
% mhmm_windowed_inference(project,outPath,trace_struct_final,bin_groups,id_string,Tres,K,w,t_window,inference_times,'savio',savio)

% Appendix 6 Figure 3 D-F: 2-state time-varying inference
w = 7;
K = 2;
id_string = 'appendix6_figure3B';
mhmm_windowed_inference(project,outPath,trace_struct_final,bin_groups,id_string,Tres,K,w,t_window,inference_times,'savio',savio)
