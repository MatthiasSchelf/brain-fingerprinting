function sensors=forwardmodel(head,sources)
% sources has dimensions [timepoints sources]
headData = load(head);
Gain = headData.Gain;
GridOrient = headData.GridOrient;
Gain=Gain(1:306,:);
constrained_gain = bst_gain_orient_standalone(Gain, GridOrient);
constrained_gain = constrained_gain(:,1:15002);
sensors = constrained_gain * sources';
end