function out = generateNoise(param)

% This function generates the noise that acts on the dynamics of the
% Fishery management problem.

% Input: param - structure containing the distributions associated with any
% source of noise in the systems

out = [random(param.v,2,1);random(param.gamma,2,1);random(param.lambda);random(param.delta,2,1)];

end
