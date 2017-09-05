function particle_pos = reDetect(im, pos, opt, param, paramSR, A_pos ,A_neg)

sz = opt.tmplsize;  
if size(im,3) == 3
  imggray = rgb2gray(im);
else
  imggray = im;
end
[~, Y, param] = affineSample(double(imggray), pos, sz, opt, param);
YY = normVector(Y);                                             
AA_pos = normVector(A_pos);
AA_neg = normVector(A_neg);
% represent each candidate with training template set
paramSR.L = length(YY(:,1));                                    
paramSR.lambda = 0.01;
beta = mexLasso(YY, [AA_pos AA_neg], paramSR);
beta = full(beta);

% the confidence value of each candidate
rec_p = sum((YY - AA_pos*beta(1:size(AA_pos,2),:)).^2);        
rec_n = sum((YY - AA_neg*beta(size(AA_pos,2)+1:end,:)).^2);
con = exp(0.01*(rec_p-rec_n)); 

[~, I] = sort(con ,'descend');
affParams = affparam2mat(param.param(:,I(1:1:opt.highest_num))); % score highest 90%

particle_pos(opt.highest_num, 2) = 0;       
particle_pos(:,1) = affParams(2,:);
particle_pos(:,2) = affParams(1,:);

end