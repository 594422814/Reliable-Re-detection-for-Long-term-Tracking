%% setup re-detector

% set thresholds for unreliability and reliability check
color_low_thres = p.color_low_thres;
color_high_thres = p.color_high_thres;
hog_low_thres = p.hog_low_thres;
hog_high_thres = p.hog_high_thres;

target = getSubwindow(im, pos, baseTargetSize);
FirstObjectMap = getColourMap(target, bg_hist, fg_hist, p.n_bins, p.grayscale_sequence);
FirstObjectMap(isnan(FirstObjectMap)) = 0;
FirstColorScore = sum(FirstObjectMap(:));
FirstColorScore_reciprocal = 1/FirstColorScore;

hogScorenum = 0;
colorScorenum = 0;
hogAver = 5;
colorAver = 1;

% set particle filter
opt = struct('numsample',50, 'affsig',[20,20,.000,.000,.000,.000]);
opt.tmplsize = [20 20];              % [height width]
sz = opt.tmplsize;                                               
opt.highest_num = round(opt.numsample * (1-0.9));    % 90% particles

params = [pos(2), pos(1) , target_sz(2) , target_sz(1) , 0.0];
param0 = [params(1), params(2), params(3)/sz(2), params(5), params(4)/params(3), 0];
p0 = params(4)/params(3); 
param0 = affparam2mat(param0);
param = [];
param.est = param0';

% set positive and negative templates
num_p = 50;                        
num_n = 200;
[A_pos, A_neg] = affineTrainG(im, sz, opt, param, num_p, num_n, p0);             
paramSR.lambda2 = 0;
paramSR.mode = 2; 
