% Script to generate TF Profiles
function main03_make_tf_input_data(varargin)

% File Paths
writePath = ['../../out/'];
% fixed data from Dubuis 2013
TFPath = '../../dat/external_data_sources/GregorData/Raw_Profiles/';
% live Bcd-GFP data from Liz and Jonathan
Bcd_path = 'E:\Jonathan\Dropbox\ZeldaLizHernan\Data\';
% Bcd_path = 'E:\Nick\Dropbox (Garcia Lab)\ZeldaLizHernan\Data\';

% ap vector
ap_vec_src = 1:1000;
ap_vec_out = 100:900;
InterpGrid = 0:20:60*60;
ap_ft = ismember(ap_vec_src,ap_vec_out);
% set max time to accept
t_lim = 50;

for i = 1:numel(varargin)
    if ischar(varargin{i})
        eval([varargin{i} ' = varargin{i+1};'])
    end
end

mkdir(writePath);
% load trace-related data
% load([dataPath '/analysis_data_mRNA.mat'])
gapFileList = dir([TFPath '*.csv']);
bcdFileList = dir([Bcd_path '*-BcdGFP-HisRFP']);

tf_input_struct = struct;
iter = 1;
for i = [1 2 4]
    % Load file and extract info
    gene_name = gapFileList(i).name;
    raw_array = csvread([TFPath gene_name],3); % load  
    % extract useful vectors
    position_vec = round(raw_array(:,2));
    sigma_vec = raw_array(:,3); % degree of cellularization
    time_vec = raw_array(:,4); % pull time column
    % sort by time (ascending)
    [time_vec, ti] = sort(time_vec); % sort   
    pt_array = raw_array(:,5:end); 
    pt_array = pt_array(ti,:);    
    sigma_vec = sigma_vec(ti);
    % filter for times within limit
    ft = time_vec<=t_lim;
    pt_array = pt_array(ft,:);
    pt_array = pt_array(:,ap_ft);
    time_vec = time_vec(ft);
    sigma_vec = sigma_vec(ft);    
    %%% Fit Approach #1: Similar to Dubuis 2013. Use Gaussian smoothing with 
    %%% sigma=5min to estimate concentration time trend 
    time_sigma = 5; % minutes 
    tf_grid_v1 = NaN(numel(InterpGrid),numel(ap_vec_out));
    % generate weight array
    wt_array = NaN(numel(sigma_vec),numel(InterpGrid));        
    for j = 1:numel(InterpGrid)
        dt_vec = time_vec - InterpGrid(j)/60;
        if min(dt_vec) > 0 || max(dt_vec) < 0
            continue
        end
        wt_vec = exp(-.5*(dt_vec/time_sigma).^2);
        wt_array = repmat(wt_vec,1,numel(ap_vec_out));
        % take weighted average of time weights
        tf_grid_v1(j,:) = nansum(wt_array.*pt_array) ./ nansum(wt_array);                        
    end    
    % subtract background
    tf_grid_v1 = tf_grid_v1 - nanmin(tf_grid_v1,[],2);    
    
    %%% Fit Approach #2: Polynomial. By fitting a second degree polynomial, I'm
    %%% imposing a prior assumption that the "true" TF levels at a given AP
    %%% position should be a smooth function in time with 1 (and only 1) peak
%     tf_grid_v2 = NaN(size(tf_grid_v1));
%     for a = 1:numel(ap_vec_out)
%         p = polyfit(time_vec, pt_array(:,a),2);
%         tf_grid_v2(:,a) = polyval(p,InterpGrid/60);
%     end
%     tf_grid_v2 = tf_grid_v2 - nanmin(tf_grid_v2,[],2); 
%     
    % Record Results 
    tf_input_struct(iter).src = gene_name;
    dashes = strfind(gene_name,'_');
    tf_input_struct(iter).TF = gene_name(dashes(1)+1:dashes(2)-1);
    tf_input_struct(iter).sigma_vec = sigma_vec;
    tf_input_struct(iter).time_vec = time_vec;
    tf_input_struct(iter).InterpGrid = InterpGrid;
    tf_input_struct(iter).ap_vec = ap_vec_out;
    tf_input_struct(iter).pt_raw = pt_array;
    tf_input_struct(iter).pt_v1 = tf_grid_v1;
%     tf_input_struct(iter).pt_v2 = tf_grid_v2;
%     
    iter = iter + 1;       
end

%%% Now generate Bcd gradient data

% data at native resolution
time_full = [];
ap_full = [];
bcd_full_norm = [];
id_full = [];
% data at interp res
time_interp = [];
ap_interp = [];
id_interp = [];
for i = 1:numel(bcdFileList)
    f_path = [Bcd_path bcdFileList(i).name];    
    if isfolder(f_path)
        load([f_path '/CompiledNuclei.mat'])
        load([f_path '/' bcdFileList(i).name '_lin.mat']);           
        RawTF = AllTracesVector;                
        % Generate AP vector
        ap_temp = [];
        time_temp = [];
        bcd_temp = []; 

        schnitz_ref = [CompiledNuclei.schnitz];
        for j = 1:numel(schnitzcells)
            cn_ind = find(schnitz_ref==j);
            if ~isempty(cn_ind)
                apv = schnitzcells(j).APpos;
                bcdv = CompiledNuclei(cn_ind).FluoMax';
                s_frames = schnitzcells(j).frames;
                cn_frames = CompiledNuclei(cn_ind).Frames;
                s_frames(apv>max(ap_vec_out)/1000|apv<min(ap_vec_out)/1000) = NaN;
                cn_frames(isnan(bcdv)) = NaN;
                % time vec
                tv = ElapsedTime(cn_frames(ismember(cn_frames,s_frames)&cn_frames>=nc14)) - ElapsedTime(nc14);
                if numel(tv) < 2
                    continue
                end          
                apv = 1000*apv(ismember(s_frames,cn_frames)&s_frames>=nc14);
                bcdv = bcdv(ismember(cn_frames,s_frames)&cn_frames>=nc14);                
                % record
                ap_temp = [ap_temp apv];
                bcd_temp = [bcd_temp bcdv];                
                time_temp = [time_temp tv];
            end
        end            
        bcd_full_norm = [bcd_full_norm bcd_temp];
        id_full = [id_full repelem(i,numel(bcd_temp))];
        time_full = [time_full time_temp];
        ap_full = [ap_full ap_temp];
    end    
end

% Generate Raw AP-Time Bcd Grid with same dims as GAP data
bcd_grid_raw = NaN(numel(time_vec),numel(ap_vec_out));
for t = 1:numel(time_vec)
    ft = round(time_full)==round(time_vec(t));
    b_vec = accumarray(round(ap_full(ft))'-min(ap_vec_out),bcd_full_norm(ft)',[numel(ap_vec_out) 1],@mean);
    b_vec(0==b_vec) = NaN;
    bcd_grid_raw(t,:) = b_vec;
end

% Because adjacent time points come from same embryos, spatial fluctuations are
% correlated in time and are not automatically smoothed out by temporal
% filtering as for the Gap Gene data. To combat this, we also apply
% gaussian filtering across space. 
ap_sigma = 5; % .5% AP
bcd_grid_smooth = NaN(size(bcd_grid_raw));
ap_mat = repmat(ap_vec_out,size(bcd_grid_raw,1),1);
for i = 1:size(bcd_grid_raw,2)
    ap_diffs = ap_mat - ap_vec_out(i);
    if ap_vec_out(i) < min(ap_full) || ap_vec_out(i) > 750 % max(ap_full) NL: anamalously high points at around 80% AP
        continue
    end
    ap_weights = exp(-.5*(ap_diffs/ap_sigma).^2);
    bcd_vec = nansum(ap_weights.*bcd_grid_raw,2) ./ nansum(ap_weights.*~isnan(bcd_grid_raw),2);
    bcd_grid_smooth(:,i) = bcd_vec;
end
% Method 1: Here I use a Gaussian kernel w/ sigma defined in time space,
% rather than invagination space
bcd_grid_v1 = NaN(numel(InterpGrid),numel(ap_vec_out));
t_mat = repmat(time_vec*60,1,size(bcd_grid_v1,2));
for i = 1:size(bcd_grid_v1,1)
    t_diffs = t_mat - InterpGrid(i);
    t_weights = exp(-.5*(t_diffs/(time_sigma*60)).^2);
    bcd_vec = nansum(t_weights.*bcd_grid_smooth) ./ nansum(t_weights.*~isnan(bcd_grid_smooth));
    bcd_grid_v1(i,:) = bcd_vec;    
end
bcd_grid_v1 = bcd_grid_v1 - nanmin(bcd_grid_v1,[],2);

tf_input_struct(iter).TF = 'Bcd';
tf_input_struct(iter).sigma_vec = NaN;
tf_input_struct(iter).time_vec = NaN;
tf_input_struct(iter).InterpGrid = InterpGrid;
tf_input_struct(iter).ap_vec = ap_vec_out;
tf_input_struct(iter).pt_raw = bcd_grid_raw;
tf_input_struct(iter).pt_v1 = bcd_grid_v1;

%%% Normalize Data. Take approach similar to one employed by Dubuis 2013
for i = 1:numel(tf_input_struct)
    % method 1
    tf1 = tf_input_struct(i).pt_v1;
%     tf1(tf1<0) = 0;
    tf1 = tf1 / (nanmax(tf1(:))-nanmin(tf1(:)));
    tf_input_struct(i).pt_v1_normed = tf1;
end

% save
save([writePath 'tf_input_struct.mat'], 'tf_input_struct')