function [responseMap, response_cf, score] = confidenceCalculate(candidate_cf, p, hann_window, hf_num, hf_den, hf_num_deep, hf_den_deep, indLayers)
% re-utilize the baseline racker

xt = getFeatureMap(candidate_cf, p.feature_type, p.cf_response_size, p.hog_cell_size);
% apply Hann window
xt_windowed = bsxfun(@times, hann_window, xt);
% compute FFT
xtf = fft2(xt_windowed);
% Correlation between filter and test patch gives the response
% Solve diagonal system per pixel.
hf = bsxfun(@rdivide, hf_num, sum(hf_den, 3) + p.lambda);
response_cf = real(ifft2(sum(conj(hf) .* xtf, 3)));
% Crop square search region (in feature pixels).
response_cf = cropFilterResponse(response_cf, floor_odd(p.norm_delta_area / p.hog_cell_size));
if p.hog_cell_size > 1 
    % Scale up to match center likelihood resolution.
    response_cf = mexResize(response_cf, p.norm_delta_area,'auto');
end

response_deep = cell(1,length(indLayers));
xtf_deep = cell(1,length(indLayers));
xt_deep  = getDeepFeatureMap(candidate_cf, hann_window, indLayers);      
for ii = 1 : length(indLayers)
   xtf_deep{ii} = fft2(xt_deep{ii});
   hf_deep = bsxfun(@rdivide, hf_num_deep{ii}, sum(hf_den_deep{ii}, 3) + p.lambda);                     
   response_deep{ii} = ensure_real( ifft2(sum( conj(hf_deep) .* xtf_deep{ii}, 3))  );
   response_deep{ii} = cropFilterResponse(response_deep{ii}, floor_odd(p.norm_delta_area / p.hog_cell_size));
   response_deep{ii} = mexResize(response_deep{ii}, p.norm_delta_area,'auto');
end

% weighting parameters [1, 0.5, 0.02] are from HCF
final_response_deep = response_deep{1} + 0.5*response_deep{2} + 0.02*response_deep{3};

responseMap = final_response_deep;
score = max(responseMap(:));

end

% We want odd regions so that the central pixel can be exact
function y = floor_odd(x)
    y = 2*floor((x-1) / 2) + 1;
end

function y = ensure_real(x)
    assert(norm(imag(x(:))) <= 1e-5 * norm(real(x(:))));
    y = real(x);
end
