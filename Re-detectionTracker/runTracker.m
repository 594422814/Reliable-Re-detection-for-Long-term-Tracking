clc;
clear all;
% RUN_TRACKER  is the external function of the tracker - does initialization and calls trackerMain
%% path
addpath('./re-detector');
addpath('./detector');                        % from Staple
addpath('./utility');

%% load video info
videoname = 'Jogging2'; 
img_path =  'sequence/Jogging2/img/';
base_path = 'sequence/';
[img_files, pos, target_sz, video_path] = load_video_info(base_path, videoname);

%% tracking-by-detection
params.grayscale_sequence = false;	       % suppose that sequence is colour
params.hog_cell_size = 4;
params.fixed_area = 150^2;                 % standard area to which we resize the target
params.n_bins = 2^5;                       % number of bins for the color histograms (bg and fg models)
params.learning_rate_pwp = 0.02;           % bg and fg color models learning rate 
params.feature_type = 'fhog';
params.inner_padding = 0.2;                % defines inner area used to sample colors from the foreground
params.output_sigma_factor = 1/16;        % standard deviation for the desired translation filter output
params.lambda = 1e-3;                      % regularization weight
params.learning_rate_cf = 0.01;            % HOG model learning rate
params.merge_factor = 0.3;                 % fixed interpolation factor - how to linearly combine the two responses
%% unreliability and reliability check
params.color_low_thres = 0.7;
params.color_high_thres = 0.8;
params.hog_low_thres = 0.6;
params.hog_high_thres = 0.7;
%% scale related
params.hog_scale_cell_size = 4;            % from DSST 
params.learning_rate_scale = 0.025;
params.scale_sigma_factor = 1/2;
params.num_scales = 33;
params.scale_model_factor = 1.0;
params.scale_step = 1.02;
params.scale_model_max_area = 32*16;
%% debugging stuff 
params.visualization = 1;                 % show output bbox on frame
params.visualization_dbg = 0;             % show also per-pixel scores, desired response and filter output

%% start trackerMain.m
im = imread([img_path img_files{1}]);

% is a grayscale sequence ?
if(size(im,3)==1)
   params.grayscale_sequence = true;
else
   params.grayscale_sequence = false;
end

params.img_files = img_files;
params.img_path = img_path;
% init_pos is the centre of the initial bounding box
params.init_pos = pos;
params.target_sz = target_sz;
[params, bg_area, fg_area, area_resize_factor] = initializeAllAreas(im, params);

if params.visualization
    params.videoPlayer = vision.VideoPlayer('Position', [100 100 [size(im,2), size(im,1)]+30]);
end

% start the actual tracking
trackerMain(params, im, bg_area, fg_area, area_resize_factor);
