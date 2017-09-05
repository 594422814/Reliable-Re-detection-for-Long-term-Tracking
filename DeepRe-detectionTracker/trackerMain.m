function [results] = trackerMain(p, im, bg_area, fg_area, area_resize_factor)
% The main framework is obtained from Staple
% Setup detector and re-detector
setupDetector;
setupRedetector;

% Main Loop
for frame = 1:num_frames
if frame>1
    
       im = imread([p.img_path p.img_files{frame}]);
       % extract patch of size bg_area and resize to norm_bg_area
       im_patch_cf = getSubwindow(im, pos, p.norm_bg_area, bg_area); 
       % compute feature map
       xt = getFeatureMap(im_patch_cf, p.feature_type, p.cf_response_size, p.hog_cell_size);
       % apply Hann window
       xt_windowed = bsxfun(@times, hann_window, xt);
       % compute FFT
       xtf = fft2(xt_windowed);
       
       % Correlation between filter and test patch gives the response
       hf = bsxfun(@rdivide, hf_num, sum(hf_den, 3) + p.lambda);
       response_cf = real(ifft2(sum(conj(hf) .* xtf, 3)));
       % Crop square search region (in feature pixels).
       response_cf = cropFilterResponse(response_cf, floor_odd(p.norm_delta_area / p.hog_cell_size));
       if p.hog_cell_size > 1
       % Scale up to match center likelihood resolution.
           response_cf = mexResize(response_cf, p.norm_delta_area,'auto');
       end
       
       xt_deep  = getDeepFeatureMap(im_patch_cf, hann_window, indLayers);      
       for ii = 1 : length(indLayers)
         xtf_deep{ii} = fft2(xt_deep{ii});
         hf_deep = bsxfun(@rdivide, hf_num_deep{ii}, sum(hf_den_deep{ii}, 3) + p.lambda);                     
         response_deep{ii} = ensure_real( ifft2(sum( conj(hf_deep) .* xtf_deep{ii}, 3))  );
         response_deep{ii} = cropFilterResponse(response_deep{ii}, floor_odd(p.norm_delta_area / p.hog_cell_size));
         response_deep{ii} = mexResize(response_deep{ii}, p.norm_delta_area,'auto');
       end
       
       % weighting parameters [1, 0.5, 0.02] are from HCF
       final_response_deep = response_deep{1} + 0.5*response_deep{2} + 0.02*response_deep{3};
       response =  final_response_deep;
              
       [row, col] = find(response == max(response(:)), 1);
       center = (1 + p.norm_delta_area) / 2;
       pos = pos + ([row, col] - center) / area_resize_factor;

     %%  calculate HOGScore and ColorScore
       hogScore = calculatePSR(response_cf);            % hogScore
       object = getSubwindow(im, pos, baseTargetSize, target_sz);     
       colorMap = getColourMap(object, bg_hist, fg_hist, p.n_bins, p.grayscale_sequence);
       colorMap(isnan(colorMap)) = 0;
       % colorScore = sum(colorMap(:))/FirstColorScore;   % colorScore
       colorScore = sum(colorMap(:))*FirstColorScore_reciprocal;   
       adaptiveUpdate;
       
      %%  "Unreliability" and "Reliability" check
       % Unreliability check 
      if  ( colorScore < color_low_thres * colorAver )||( hogScore < hog_low_thres * hogAver )     
           disp( [num2str(frame),'th Frame. Searching if there exists a better choice.']); 
            % coarse localization
           particle_pos = reDetect(im, pos, opt, param, paramSR, A_pos ,A_neg);          
           responseFinal = cell(1,opt.highest_num);  
           % discard 90% of the particles
           responseCF = cell(1,opt.highest_num); 
           score(1,opt.highest_num) = 0;
           dist(1,opt.highest_num) = 0;
           d(1,opt.highest_num) = 0;
           candidate_pos(opt.highest_num, 2) = 0;
           % distance penalize (in case of target too small)
           dist_param = max(sum(target_sz)/2, 35); 
           
          for i = 1 : opt.highest_num
              candidate_cf = getSubwindow( im, particle_pos(i,:), p.norm_bg_area, bg_area);
              [responseFinal{i}, responseCF{i}, score(i)] = confidenceCalculate(candidate_cf, p ,hann_window, hf_num, hf_den, hf_num_deep, hf_den_deep, indLayers);
             
              [row, col] = find( responseFinal{i} == max(responseFinal{i}(:)), 1);    
              candidate_pos(i,:) = particle_pos(i,:) + ([row, col] - center)/area_resize_factor;
              % distance penalize
              dist(i) = sqrt( (candidate_pos(i,1) - pos(1))^2 + (candidate_pos(i,2) - pos(2))^2);
              d(i) = cosd( (10 * dist(i))/dist_param ); 
              score(i) = score(i)*d(i);
          end  
           
          % Choose the best candidate
          [~,ID] = max(score);                                              
          % Calculate candidate hogScore and colorScore
          CanHOGS = calculatePSR(responseCF{ID});
          object = getSubwindow(im, candidate_pos(ID,:), baseTargetSize, target_sz);
          objectScore = getColourMap(object, bg_hist, fg_hist, p.n_bins, p.grayscale_sequence);
          objectScore(isnan( objectScore )) = 0;
          CanColorS = sum(objectScore(:))*FirstColorScore_reciprocal; 
          
          % "Reliability" check
          if (d(ID)*CanColorS > color_high_thres*colorAver)&&(d(ID)*CanHOGS > hog_high_thres*hogAver)   
              pos = candidate_pos(ID,:);
          end  
      end
        rect_position = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])]; 
          
      %% SCALE SPACE SEARCH
        im_patch_scale = getScaleSubwindow(im, pos, base_target_sz, scale_factor * scale_factors, scale_window, scale_model_sz, p.hog_scale_cell_size, hann_window);
        xsf = fft(im_patch_scale,[],2);
        scale_response = real(ifft(sum(sf_num .* xsf, 1) ./ (sf_den + p.lambda) ));
        recovered_scale = ind2sub(size(scale_response),find(scale_response == max(scale_response(:)), 1));
        %set the scale
        scale_factor = scale_factor * scale_factors(recovered_scale);

        if scale_factor < min_scale_factor
            scale_factor = min_scale_factor;
        elseif scale_factor > max_scale_factor
            scale_factor = max_scale_factor;
        end
               
        % use new scale to update bboxes for target, filter, bg and fg models
        target_sz = round(base_target_sz * scale_factor);
        p.avg_dim = sum(target_sz)/2;
        bg_area = round(target_sz + p.avg_dim * p.padding);   
        if(bg_area(2)>size(im,2)),  bg_area(2)=size(im,2)-1;    end
        if(bg_area(1)>size(im,1)),  bg_area(1)=size(im,1)-1;    end

        bg_area = bg_area - mod(bg_area - target_sz, 2);
        fg_area = round(target_sz - p.avg_dim * p.inner_padding);
        fg_area = fg_area + mod(bg_area - fg_area, 2);
        % Compute the rectangle with (or close to) params.fixed_area and same aspect ratio as the target bboxgetScaleSubwindow
        area_resize_factor = sqrt(p.fixed_area/prod(bg_area));

        if p.visualization_dbg==1
            mySubplot(2,1,5,1,im_patch_cf,'FG+BG','gray');
            mySubplot(2,1,5,2,likelihood_map,'obj.likelihood','parula');
            mySubplot(2,1,5,3,response_cf,'CF response','parula');
            mySubplot(2,1,5,4,response_pwp,'center likelihood','parula');
            mySubplot(2,1,5,5,response,'merged response','parula');
            drawnow
        end

end
  %% TRAINING
    % extract patch of size bg_area and resize to norm_bg_area
    im_patch_bg = getSubwindow(im, pos, p.norm_bg_area, bg_area);
    % compute feature map, of cf_response_size
    xt = getFeatureMap(im_patch_bg, p.feature_type, p.cf_response_size, p.hog_cell_size);
    % apply Hann window
    xt = bsxfun(@times, hann_window, xt);
    % compute FFT
    xtf = fft2(xt);
    % FILTER UPDATE
    % Compute expectations over circular shifts, therefore divide by number of pixels.
    new_hf_num = bsxfun(@times, conj(yf), xtf) / prod(p.cf_response_size);
    new_hf_den = (conj(xtf) .* xtf) / prod(p.cf_response_size);

    xt_deep  = getDeepFeatureMap(im_patch_bg, hann_window, indLayers);
    for  ii = 1 : numLayers
       xtf_deep{ii} = fft2(xt_deep{ii});
       new_hf_num_deep{ii} = bsxfun(@times, conj(yf), xtf_deep{ii}) / prod(p.cf_response_size);
       new_hf_den_deep{ii} = (conj(xtf_deep{ii}) .* xtf_deep{ii}) / prod(p.cf_response_size);
    end
    
    if frame == 1
        % first frame, train with a single image
        hf_den = new_hf_den;
        hf_num = new_hf_num;
        for ii=1 : numLayers 
           hf_den_deep{ii} = new_hf_den_deep{ii};
           hf_num_deep{ii} = new_hf_num_deep{ii};
        end
    else
        % subsequent frames, update the model by linear interpolation
        hf_den = (1 - p.learning_rate_cf) * hf_den + p.learning_rate_cf * new_hf_den;
        hf_num = (1 - p.learning_rate_cf) * hf_num + p.learning_rate_cf * new_hf_num;
        for ii= 1 : numLayers
           hf_den_deep{ii} = (1 - p.learning_rate_cf) * hf_den_deep{ii} + p.learning_rate_cf * new_hf_den_deep{ii};
           hf_num_deep{ii} = (1 - p.learning_rate_cf) * hf_num_deep{ii} + p.learning_rate_cf * new_hf_num_deep{ii};
        end
        % BG/FG MODEL UPDATE   patch of the target + padding
        [bg_hist, fg_hist] = updateHistModel(new_pwp_model, im_patch_bg, bg_area, fg_area, target_sz, p.norm_bg_area, p.n_bins, p.grayscale_sequence, bg_hist, fg_hist, p.learning_rate_pwp);
        % update positive and negative templates every 5 frames 
        if (rem(frame,5) == 0)&&(colorScore > color_high_thres * colorAver)&&(hogScore > hog_high_thres * hogAver)
            params = [pos(2), pos(1) , target_sz(2) , target_sz(1) , 0.0];
            param0 = [params(1), params(2), params(3)/sz(2), params(5), params(4)/params(3), 0];
            p0 = params(4)/params(3); 
            param0 = affparam2mat(param0);
            param.est = param0';
            [A_pos] = updatePos(im, sz, param, num_p, A_pos);
            [A_neg] = updateNeg(im, sz, param, num_n, p0);
        end
    end
    
   %% SCALE UPDATE
    im_patch_scale = getScaleSubwindow(im, pos, base_target_sz, scale_factor*scale_factors, scale_window, scale_model_sz, p.hog_scale_cell_size, hann_window);
    xsf = fft(im_patch_scale,[],2);
    new_sf_num = bsxfun(@times, ysf, conj(xsf));
    new_sf_den = sum(xsf .* conj(xsf), 1);
    if frame == 1,
        sf_den = new_sf_den;
        sf_num = new_sf_num;
    else
        sf_den = (1 - p.learning_rate_scale) * sf_den + p.learning_rate_scale * new_sf_den;
        sf_num = (1 - p.learning_rate_scale) * sf_num + p.learning_rate_scale * new_sf_num;
    end
    
    % update bbox position
    if frame==1, rect_position = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])]; end  
    rect_position_padded = [pos([2,1]) - bg_area([2,1])/2, bg_area([2,1])];  

   %% VISUALIZATION
    if p.visualization == 1
        if isToolboxAvailable('Computer Vision System Toolbox')
            im = insertShape(im, 'Rectangle', rect_position, 'LineWidth', 4, 'Color', 'red');
            im = insertShape(im, 'Rectangle', rect_position_padded, 'LineWidth', 4, 'Color', 'yellow');
            % Display the annotated video frame using the video player object.
            step(p.videoPlayer, im);
       else
            figure(1)
            imshow(im)
            rectangle('Position',rect_position, 'LineWidth',2, 'EdgeColor','r');
            rectangle('Position',rect_position_padded, 'LineWidth',2, 'LineStyle','-', 'EdgeColor','y');
            drawnow
       end
    end
end

end


% We want odd regions so that the central pixel can be exact
function y = floor_odd(x)
    y = 2*floor((x-1) / 2) + 1;
end

function y = ensure_real(x)
    assert(norm(imag(x(:))) <= 1e-5 * norm(real(x(:))));
    y = real(x);
end

function hogScore = calculatePSR(response_cf)
    cf_max = max(response_cf(:));
    cf_average = mean(response_cf(:));
    cf_sigma = sqrt(var(response_cf(:)));
    hogScore = (cf_max - cf_average)/cf_sigma;
end