function param = ComputeTransition(TypeVectorField,NumberOfPartitions,NumberOfMonteCarlo)

switch TypeVectorField
    case 'TCL'
        %% Model parameters

        R = 2; % Thermal resistence in Celsius/kW
        C = 2; % Thermal capacity ini kW/Celsius
        P = 14; % Range of energy transfer to or from the thermal mass in kW
        eta = 0.7; % Control efficiency
        h = 5/60; % Discretization time in hours

        alpha = exp(-h/(C*R));

        Al = 19; % Lower bound on the safe set
        Ah = 22; % Upper bound on the safe set

        mu = 0; sigma = 100*0.25^2; % Estimate on the mean and variance
        W = [-0.5*sqrt(sigma/12),0.5*sqrt(sigma/12)];  % Support of the distribution

        theta = 32; % Environment temperature
        p = 0.05; % success probability

        param.alpha = alpha;
        param.theta = theta;
        param.eta = eta;
        param.R = R;
        param.P = P;
        param.C = C;
        param.h = h;
        param.p = p;
        param.SafeSet = [Al,Ah];

        pd = makedist('Normal','mu',mu,'sigma',sigma);
        pd_truncated = truncate(pd,W(1),W(2));

        param.w = pd_truncated;

        param.NumberOfPartitions = NumberOfPartitions;
        InputPartition = generateInputPartition([],'TCL'); % Generate a vector with all possible combinations of inputs


        Grid = StatePartition(NumberOfPartitions,param.SafeSet,'TCL'); % Generate the partition of the state space
        [List,ListX] = Grid.createList(InputPartition); % List containing labels for the discrete states of the discretazation

        param.List = List;
        param.ListX = ListX;
        param.sizePartition = Grid.getSizePartition;

        param.MC = NumberOfMonteCarlo;
        
        TransitionProb = EstimateTransition(Grid,InputPartition,param); 

        % Generating the transition probability
        param.TransitionProb = TransitionProb;

    otherwise
        NotImplemented();
end




end