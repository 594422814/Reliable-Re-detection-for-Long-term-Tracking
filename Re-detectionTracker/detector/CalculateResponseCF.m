function [ response_cf ] = CalculateResponseCF(im,pos,p,bg_area,hann_window,model_xf,model_alphaf)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % extract patch of size bg_area and resize to norm_bg_area
    im_patch_cf = getSubwindow(im, pos, p.norm_bg_area, bg_area);
    switch CFtype
        case 'HOGCN'
           xf = get_features(im_patch_cf, p.feature_type, p.cf_response_size, p.hog_cell_size, w2c);  
        case 'HOG'
           xf = getFeatureMap(im_patch_cf, p.feature_type, p.cf_response_size, p.hog_cell_size);
    end
    % apply Hann window
    xf_windowed = bsxfun(@times, hann_window, xf);
    % compute FFT
    xf = fft2(xf_windowed);
    kxf = sum(xf .* conj(model_xf), 3) / numel(xf);
    response_cf = real(ifft2(kxf .* model_alphaf /p.lambda1 ));
    % Crop square search region (in feature pixels).
    response_cf = cropFilterResponse(response_cf, floor_odd(p.norm_delta_area / p.hog_cell_size));
    if p.hog_cell_size > 1
        % Scale up to match center likelihood resolution.
        response_cf = mexResize(response_cf, p.norm_delta_area,'auto');
    end
    response_cf = response_cf./max(response_cf(:));
            
end

% We want odd regions so that the central pixel can be exact
function y = floor_odd(x)
    y = 2*floor((x-1) / 2) + 1;
end
