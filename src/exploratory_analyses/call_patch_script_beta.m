% Script to call make_patch_movies function for mean rate, frac on, and
% duration
clear
close all
% set file paths
PreProcPath = 'E:\Eve2MSE\Data\PreProcessedData\';
PostProcPath = 'E:\Eve2MSE\Data\DynamicsResults/';
%Prefix = '2019-03-22-eve2-MSE-Hb_MCP2f_02';
nc = 14;

Prefix_index = {'2019-03-22-eve2-MSE-Hb_MCP2f_02','2019-03-27-eve2-MSE-Gt_MCP2f_01','2019-03-02-eve2-MSE-Wt_MCP2f_01'};
color_index = {[234 194 100],[11 190 190],[85 169 116]};

for i=2%:numel(Prefix_index)
    Prefix = Prefix_index{i};
    color = color_index{i};
    make_patch_movies_beta(Prefix, PreProcPath, PostProcPath, nc, 'instantRate','visible','last_time',50, 'set_color', color)
    %make_patch_movies_beta(Prefix, PreProcPath, PostProcPath, nc)
end
% call function for fraction on
%make_patch_movies_beta(Prefix, PreProcPath, PostProcPath, nc, 'instantRate','visible','last_time',40)

% call function for mean rate
%make_patch_movies_beta(Prefix, PreProcPath, PostProcPath, nc)

