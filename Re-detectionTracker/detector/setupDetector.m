%% set up detector (from Staple)

pos = p.init_pos;
target_sz = p.target_sz;
baseTargetSize = target_sz;
num_frames = numel(p.img_files);

% patch of the target + padding
patch_padded = getSubwindow(im, pos, p.norm_bg_area, bg_area);

% initialize hist model
new_pwp_model = true;
[bg_hist, fg_hist] = updateHistModel(new_pwp_model, patch_padded, bg_area, fg_area, target_sz, p.norm_bg_area, p.n_bins, p.grayscale_sequence);
new_pwp_model = false;
% Hann (cosine) window
 hann_window = single(hann(p.cf_response_size(1)) * hann(p.cf_response_size(2))');
% gaussian-shaped desired response, centred in (1,1)
% bandwidth proportional to target size
output_sigma = sqrt(prod(p.norm_target_sz)) * p.output_sigma_factor / p.hog_cell_size;
y = gaussianResponse(p.cf_response_size, output_sigma);
yf = fft2(y);

%% SCALE ADAPTATION INITIALIZATION
% Code from DSST
scale_factor = 1;
base_target_sz = target_sz;
scale_sigma = sqrt(p.num_scales) * p.scale_sigma_factor;
ss = (1:p.num_scales) - ceil(p.num_scales/2);
ys = exp(-0.5 * (ss.^2) / scale_sigma^2);
ysf = single(fft(ys));
if mod(p.num_scales,2) == 0
    scale_window = single(hann(p.num_scales+1));
    scale_window = scale_window(2:end);
else
    scale_window = single(hann(p.num_scales));
end;

ss = 1:p.num_scales;
scale_factors = p.scale_step.^(ceil(p.num_scales/2) - ss);

if p.scale_model_factor^2 * prod(p.norm_target_sz) > p.scale_model_max_area
    p.scale_model_factor = sqrt(p.scale_model_max_area/prod(p.norm_target_sz));
end

scale_model_sz = floor(p.norm_target_sz * p.scale_model_factor);
% find maximum and minimum scales
min_scale_factor = p.scale_step ^ ceil(log(max(5 ./ bg_area)) / log(p.scale_step));
max_scale_factor = p.scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ target_sz)) / log(p.scale_step));
