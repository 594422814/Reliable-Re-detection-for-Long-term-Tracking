function [wimgs, Y, param] = affineSample(frm, pos, sz, opt, param)
% draw N candidates with particle filter

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
