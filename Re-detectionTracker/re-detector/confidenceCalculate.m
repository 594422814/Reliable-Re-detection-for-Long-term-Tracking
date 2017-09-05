function [responseMap, response_cf, score] = confidenceCalculate(candidate_pwp, candidate_cf, bg_hist, fg_hist, p, hann_window, hf_num, hf_den)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
colormap = getColourMap(candidate_pwp , bg_hist, fg_hist, p.n_bins, p.grayscale_sequence);
colormap(isnan(colormap)) = 0;
response_pwp = getCenterLikelihood(colormap, p.norm_target_sz);

xt = getFeatureMap(candidate_cf, p.feature_type, p.cf_response_size, p.hog_cell_size);
% apply Hann window
xt_windowed = bsxfun(@times, hann_window, xt);
% compute FFT
xtf = fft2(xt_windowed);
% Correlation between filter and test patch gives the response
% Solve diagonal system per pixel.
hf = bsxfun(@rdivide, hf_num, sum(hf_den, 3)+p.lambda);
response_cf = real(ifft2(sum(conj(hf) .* xtf, 3)));
% Crop square search region (in feature pixels).
response_cf = cropFilterResponse(response_cf, floor_odd(p.norm_delta_area / p.hog_cell_size));
if p.hog_cell_size > 1 
    % Scale up to match center likelihood resolution.
    response_cf = mexResize(response_cf, p.norm_delta_area,'auto');
end

responseMap = (1-p.merge_factor) * response_cf + p.merge_factor * response_pwp;
score = max(responseMap(:));

end

% We want odd regions so that the central pixel can be exact
function y = floor_odd(x)
    y = 2*floor((x-1) / 2) + 1;
end
