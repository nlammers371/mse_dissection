% Master script to call mHMM inference for various conditions of interest
clear
close all
% First define variables that are common across all inference runs
project = 'mse_comparison_lateralML'; % Project identifier
Tres = 20; % Time Resolution
t_window = 60;
inference_times = 30;
savio = 0;
if savio
    outPath = '/global/scratch/nlammers/mse_inference_results/'; %hmmm_data/inference_out/';
else    
    outPath = '../../out/ss_inference_results/';
end
mkdir(outPath);
% load inference traces 
dataPath = ['../../dat/' project '/']; %Path to raw data
dataName = ['qc_nucleus_struct.mat'];
load([dataPath dataName]);
nucleus_struct = nucleus_struct([nucleus_struct.qc]==1);
for i = 1:numel(nucleus_struct)
    nucleus_struct(i).ap_id = round(100*mean(nucleus_struct(i).ap_vector_interp));
    nucleus_struct(i).alpha_frac = 1302/6444;
end

% designate inference regions
% ap_ref_cell = {-5:-4,-3,-2,-1,0,1,2,3,4:5}; 
ap_ref_cell = {38:39};%26:29,30:33,34:37,38:49,40:41,42:43,44:47,48:52};%{20:25,26:30,31:32,
% ap_ref_cell = {37:38,39:40};
gtype_id_index = 3;
bin_groups = ap_ref_cell;%{-3};%,-2,-1,0,1,2,3};%:1,2:3,4:5};
if savio
    %Get environment variable from job script
    savio_groups = {str2num(getenv('SLURM_ARRAY_TASK_ID'))};    
    bin_groups = cell(1,length(savio_groups));
    for i =1:length(bin_groups)
        bin_groups{i} = ap_ref_cell{savio_groups{i}};
    end
end  
% 
% % Figure 5: 3 State time-averaged inference 
w = 6;
K = 3;
n_boots = 5;
% id_string = {'Hb'};
id_string = {'Hb'};
mhmm_windowed_inference(project,outPath,nucleus_struct,bin_groups,id_string,...
    Tres,K,w,t_window,inference_times,gtype_id_index,'savio',savio,'n_bootstrap',n_boots)



