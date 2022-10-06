function out = generateNoise(param,type)

% This function generates the noise that acts on the dynamics of the
% Fishery management problem and TCL problem in Insoon Yang's paper.

% Input: param - structure containing the distributions associated with any
% source of noise in the systems
%        type -- this is the type of noise we want to generate, depending on the vector field 
%        in use.
%

switch type
    case 'Fishery'
        out = [random(param.v,2,1);random(param.gamma,2,1);random(param.lambda);random(param.delta,2,1)];
    case 'TCL'
        out = random(param.w);
    case 'ChainInt'
        out = random(param.w,2,1);
    case 'CarPole'
        out = random(param.w,4,1);
    case 'CarPoleNL'
        out = random(param.w,4,1);
    otherwise
        NotImplemented();
        
end
