function [A_pos] = updatePos(img, sz, param, num_p, A_pos)

if size(img,3)==3
    img	= double(rgb2gray(img));
else
    img	= double(img);
end

n = num_p;                     % Sampling Number
param.param0 = zeros(6,n);     % Affine Parameter Sampling
param.param = zeros(6,n);
param.param0 = repmat(affparam2geom(param.est(:)), [1,n]);
randMatrix = randn(6,n);
sigma = [1, 1, .000, .000, .000, .000];
param.param = param.param0 + randMatrix.*repmat(sigma(:),[1,n]);

o = affparam2mat(param.param);     % Extract or Warp Samples which are related to above affine parameters
wimgs = warpimg(img, o, sz);

m = prod(sz);

for i = 1: n/2
    A_pos(:,i) = reshape(wimgs(:,:,i), m, 1);
end

end

