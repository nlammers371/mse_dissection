% Master script to call all data cleaning functions
%run from inside the data_pipeline folder
clear 
close all

%%%%%%%%%%%%%%%%%% Call function to compile traces
% Set ID variables
% folderPath = 'E:/Nick/Dropbox (Garcia Lab)/mHMM/weka/';

project = 'mse_comparison_lateralML';
folderPath = 'E:\Eve2MSE\Data\DynamicsResults\';
keyword = 'MCP2f';
ap_shift = 37;
% keyword = '20sec'; % Keyword to ensure only sets from current project are pulled
% include_vec = [9:13 19 20 22:24 26]; %data set numbers to include
% Tres_interp = 20; % interpolated time resolution
% n_boots = 100;
% boot_arg = 'n_boots';
% specify which scripts to run
script_tape = [1,2,3,4];


function_cell = {'main01_compile_data(project,folderPath,keyword)',...
                 'main02_quality_control_beta(project)',...
                 'main03_make_mRNA_analaysis_set(project)',...
                 'main04_join_trace_and_tf_data(project,ap_shift)'};

for s = script_tape    
    disp(['Running script number ' num2str(s) '...'])     
    eval(function_cell{s})           
end
