function [wimgs, Y, param] = affineSample(frm, pos, sz, opt, param)
% function [wimgs Y param] = affineSample(frm, sz, opt, param)
% draw N candidates with particle filter
% input --- 
% frm: the image of the current frame
% sz: the size of the tracking window
% opt: initial parameters
% param: the affine parameters
% output ---
% wimgs: the N candidate images (matrix)
% Y: the N candidate images (vector)
% param: the affine parameters
%*************************************************************
%% Copyright (C) Wei Zhong.
%% All rights reserved.
%% Date: 05/2012

n = opt.numsample;                  % Sampling Number
param.param0 = zeros(6,n);          % Affine Parameter Sampling
param.param = zeros(6,n);
param.param0 = repmat(affparam2geom(param.est(:)), [1,n]);
randMatrix = randn(6,n);
param.param = param.param0 + randMatrix.*repmat(opt.affsig(:),[1,n]);
o = affparam2mat(param.param);    % Extract or Warp Samples which are related to above affine parameters
% set one particle at the the original position
o(1,n) = pos(2);
o(2,n) = pos(1);

wimgs = warpimg(frm, o, sz);

m = prod(opt.tmplsize);             % vectorization
Y = zeros(m, n);
for i = 1:n
    Y(:,i) = reshape(wimgs(:,:,i), m, 1);
end
