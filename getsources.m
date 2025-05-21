function sources=getsources(head,data,atlas,gainonly)
% PLEASE NOTE THAT DATA IS RESAMPLED HERE. COMMENT IF UNWANTED
if nargin<4
    gainonly=0;
end
headData = load(head);
Gain = headData.Gain;
GridOrient = headData.GridOrient;
Gain=Gain(1:306,:);
constrained_gain = bst_gain_orient_standalone(Gain, GridOrient);
constrained_gain = constrained_gain(:,1:15002);
if ischar(data)
    dataData = load(data);
    F = dataData.F(1:306,201:end);
    F = resample(F', 1, 3)';
elseif isstruct(data)
    F=data.F;
end
npoints=size(F,2);
if gainonly
    sources = (constrained_gain' * F)';
else
    [U, S, V] = svd(constrained_gain, 'econ');
    s = diag(S);
    s_filt = s ./ (s.^2 + 0.05);
    K = V * diag(s_filt) * U';  % MNE kernel
    % Compute variance of each source (diagonal of KK')
    kernel_power = sqrt(sum(K.^2, 2));  % same as sqrt(diag(K*K'))
    % Normalize each row to unit variance (sLORETA step)
    kernel_sLORETA = bsxfun(@rdivide, K, kernel_power + eps);  % eps for numerical stability
    sources = (kernel_sLORETA * F)';
end
if ~isempty(atlas)
    nregions=length(atlas);
    sources_avg=zeros(npoints,nregions);
    for i=1:nregions
        sources_avg(:,i)=mean(sources(:,atlas{i}),2);
    end
    sources=sources_avg;
end

end