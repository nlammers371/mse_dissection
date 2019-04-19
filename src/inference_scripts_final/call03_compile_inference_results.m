% Script to call inference result compiler 
clear
close all

% id variablepl
project = 'mse_comparison_lateralML';

% Figure 5: 3 State time-averaged inference 
w = 6;
K = 3;
figure_id_cell = {'Hb','Wt','Gt'};
inference_type = 'ss_inference_results';

for i = 1:numel(figure_id_cell)    
    compile_inference_results(project,inference_type,figure_id_cell{i},K,w)
end