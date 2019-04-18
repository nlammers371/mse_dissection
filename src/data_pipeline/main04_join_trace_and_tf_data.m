% Short script to combine tf data with trace data

function main04_join_trace_and_tf_data(project,ap_shift,varargin)
dataPath = ['../../dat/' project '/'];
% ap_shift = 37; % offset between data sets (hard code for now)
load('../../out/tf_input_struct.mat')
load([dataPath 'qc_nucleus_struct.mat'])

ap_vec = tf_input_struct(1).ap_vec - ap_shift;
time_vec = tf_input_struct(1).InterpGrid;

% add tf level array for each time step for each nucleus
for i = 1:numel(nucleus_struct)
    nc_time = nucleus_struct(i).time_interp;
    nc_ap = round(nucleus_struct(i).ap_vector_interp*1000);
    tf_array = NaN(numel(nc_time),numel(tf_input_struct));
    for j = 1:numel(tf_input_struct)
        tf_grid = tf_input_struct(j).pt_v1_normed;
        %AM: Edited this section
        [~,col] = ismember(nc_ap,ap_vec);
        [~,row] = ismember(nc_time,time_vec);
        indices = sub2ind(size(tf_grid),row,col);
        tf_array(:,j) = tf_grid(indices);
    end
    nucleus_struct(i).tf_array = tf_array;
end
save([dataPath 'nucleus_struct_protein.mat'],'nucleus_struct')