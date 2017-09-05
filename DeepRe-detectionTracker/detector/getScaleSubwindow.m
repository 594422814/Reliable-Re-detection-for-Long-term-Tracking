function out = getScaleSubwindow(im, pos, base_target_sz, scale_factors, scale_window, scale_model_sz, hog_scale_cell_size, cos_window)

% code from DSST

num_scales = length(scale_factors);

for s = 1:num_scales
    patch_sz = floor(base_target_sz * scale_factors(s));
    %make sure the size is not to small
    patch_sz = max(patch_sz, 2);
 
    xs = floor(pos(2)) + (1:patch_sz(2)) - floor(patch_sz(2)/2);
    ys = floor(pos(1)) + (1:patch_sz(1)) - floor(patch_sz(1)/2);
    
    %check for out-of-bounds coordinates, and set them to the values at
    %the borders
    xs(xs < 1) = 1;
    ys(ys < 1) = 1;
    xs(xs > size(im,2)) = size(im,2);
    ys(ys > size(im,1)) = size(im,1);
    
    %extract image
    im_patch = im(ys, xs, :);
    
    % resize image to model size
    im_patch_resized = mexResize(im_patch, scale_model_sz, 'auto');
    
%     feat = getScaleFeature(im_patch_resized, cos_window );
%     temp = feat;
    % extract scale features
    temp_hog = fhog(single(im_patch_resized), hog_scale_cell_size);
    temp = temp_hog(:,:,1:31);
    
%     	x = double(fhog(single(im_patch_resized) / 255, hog_scale_cell_size,9));
% 		x(:,:,end) = [];  %remove all-zeros channel ("truncation feature")
% 		sz = size(x);
% 		im_patch = imresize(im_patch_resized, [sz(1) sz(2)]);
% 		out_npca = get_feature_map(im_patch, 'gray', w2c);
% 		out_pca = get_feature_map(im_patch, 'cn', w2c);
% % 		out_pca = reshape(temp_pca, [prod(sz), size(temp_pca, 3)]);
% 		x = cat(3,x,out_npca);
% 		x=cat(3,x,out_pca);
%         temp = x(:,:,1:41);

        
    if s == 1
        out = zeros(numel(temp), num_scales, 'single');
    end
    
    % window
    out(:,s) = temp(:) * scale_window(s);
end
