% Master script to call input/output analysis scripts
clear
close all

project = 'mHMMeve2_weka_inf_2018_05_07';
% build control variables
call_build_functions = 0;
build_script_tape = [1,2,3,4]; % scripts to call
% inference control variables
call_inference_functions = 1;
inference_id_cell = {'main_text','appendix'};
inference_script_tape = [1,2,3]; % scripts to call
boot_arg = 'n_boots';
n_boots = 100;

if call_build_functions
    % call data-building functions   
    build_script_cell = {'main01_make_tf_input_data(project)',...
                         'main02_make_and_align_stripe2_profile(project)',...
                         'main03_compile_regression_set(project)',...
                         'main04_make_regression_weights(project)'};

    for s = 1:numel(build_script_tape)
        disp(['Running build script ' num2str(build_script_tape(s)) '...'])
        eval(build_script_cell{build_script_tape(s)})
        disp('Done.')
    end
end

if call_inference_functions
    % call inference functions    
    inference_script_cell = {'[use_weights,constrained_inference,ap_raw_flag] = main05_run_regressions(project,inference_id,boot_arg,n_boots)',...
                             'main06_add_trace_model_motifs(project,inference_id,use_weights,constrained_inference,ap_raw_flag)',...
                             'main07_make_summary_figs(project,inference_id,use_weights,constrained_inference,ap_raw_flag)'};
    
    for i = 1:numel(inference_id_cell)
        tic
        inference_id = inference_id_cell{i};
        for j = 1:numel(inference_script_tape)
            disp(['Running inference script ' num2str(inference_script_tape(j)) '...'])
            eval(inference_script_cell{inference_script_tape(j)})
            disp('Done.')
        end
        toc
    end
end