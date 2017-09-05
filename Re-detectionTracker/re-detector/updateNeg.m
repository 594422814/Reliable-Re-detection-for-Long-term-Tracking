function [X_neg] = updateNeg(img, sz, param, num_n, p0)

%*************************************************************
%% Copyright (C) Wei Zhong.
%% All rights reserved.
%% Date: 05/2012

if size(img,3)==3
    img	= double(rgb2gray(img));
else
    img	= double(img);
end

n = num_n;                                                   % Sampling Number
param.param0 = zeros(6,n);                                   % Affine Parameter Sampling
param.param = zeros(6,n);
param.param0 = repmat(affparam2geom(param.est(:)), [1,n]);
randMatrix = randn(6,n);
sigma = [round(sz(2)*param.est(3)), round(sz(1)*param.est(3)*p0), .000, .000, .000, .000];
param.param = param.param0 + randMatrix.*repmat(sigma(:),[1,n]);

back = round(sigma(1)/5);
center = param.param0(1,1);
left = center - back;
right = center + back;
nono = param.param(1,:)<=right&param.param(1,:)>=center;
param.param(1,nono) = right;
nono = param.param(1,:)>=left&param.param(1,:)<center;
param.param(1,nono) = left;

back = round(sigma(2)/5);
center = param.param0(2,1);
top = center - back;
bottom = center + back;
nono = param.param(2,:)<=bottom&param.param(2,:)>=center;
param.param(2,nono) = bottom;
nono = param.param(2,:)>=top&param.param(2,:)<center;
param.param(2,nono) = top;

o = affparam2mat(param.param);     % Extract or Warp Samples which are related to above affine parameters
wimgs = warpimg(img, o, sz);

m = prod(sz);
X_neg = zeros(m, n);
for i = 1: n
    X_neg(:,i) = reshape(wimgs(:,:,i), m, 1);
end

end
