function [] = PlotResults(OuputOfMonteCarloSimulationFunc)

NumberOfPoints = OuputOfMonteCarloSimulationFunc{end}.PartitionSize;
m = OuputOfMonteCarloSimulationFunc{end}.m;
ep = OuputOfMonteCarloSimulationFunc{end}.AmbiguityParam.ep;
rhoMu = OuputOfMonteCarloSimulationFunc{end}.AmbiguityParam.rhoMu;
rhoSigma = OuputOfMonteCarloSimulationFunc{end}.AmbiguityParam.rhoSigma;

TypeOfVectorField = OuputOfMonteCarloSimulationFunc{1}.TypeOfVectorField;
Grid = OuputOfMonteCarloSimulationFunc{end}.Partition;
AmbiguityParam = OuputOfMonteCarloSimulationFunc{end}.AmbiguityParam;

L = length(OuputOfMonteCarloSimulationFunc) - 1;

switch TypeOfVectorField
    case 'TCL'
        X = Grid.getValues.Partition.X;
        figure;
        hold on
        Name = cell(L,1);
        for i =1:L
            grid on;
            box on;
            Name{i} = OuputOfMonteCarloSimulationFunc{i}.AmbiguityType;
            ValueFunction = OuputOfMonteCarloSimulationFunc{i}.ValueFunction;
            plot(X,ValueFunction,'LineWidth',1.5);
            %plot(X,eval(Name{i}).ValueFunction(1:end-1,end-3),'LineWidth',1.5);
        end
        legend(Name)
        
        TitleString = sprintf('TCL with Partition = %d, m=%d, ep=%.2f, rhoMu = %.2f, rhoSigma = %.2f',NumberOfPoints,m,ep,rhoMu,rhoSigma);
        title(TitleString)
        
        figure;
        hold on
        Name = cell(L+1,1);
        for i =1:L
            grid on;
            box on;
            Name{i} = OuputOfMonteCarloSimulationFunc{i}.AmbiguityType;
            ValueFunction = OuputOfMonteCarloSimulationFunc{i}.ValueFunction;
            EmpValueFunction = OuputOfMonteCarloSimulationFunc{i}.EmpiricalValueFunc;
            plot(ValueFunction,EmpValueFunction,'o','LineWidth',1.5);
            %plot(X,eval(Name{i}).ValueFunction(1:end-1,end-3),'LineWidth',1.5);
        end
        Name{L+1} = '';
        plot(0:0.05:1,0:0.05:1,'--','LineWidth',1.5)
        legend(Name)
        title(TitleString)
        
        
    case 'Fis'
        Partition = Grid.getValues.Partition;
        Name = cell(L,1);
        tempVariable(:) = whos('-regexp','ValueFunc*');
        
        h = figure;
        hold on;
                
        InitialState_X3 = 750;
        InitialState_X3 = Partition.X3(1,1,find(Partition.X3(1,1,:) <= InitialState_X3,1,'last'));
        
        
        for i=1:L
            Name{i} = tempVariable(1,i).name;
            [XX,YY,ProjValueFunc,OptPolicy] = PlotFishery(InitialState_X3,eval(Name{i}).ValueFunction(:,1),eval(Name{i}).OptInput(:,1),Partition);
            levels = [0.1,0.2,0.4,0.5,0.6,0.8,0.85,0.9,1];
            figure
            contour(XX,YY,ProjValueFunc,levels,'ShowText','on','LineWidth',2);           
        end
        
    otherwise
        NotImplemented();
end

end
% Functions created to assist the plotting of the obtained results
function [X,Y,ProjValueFunction,OptPolicy] = PlotFishery(Height,ValueFunction,OptInput,...
    Partition)

Indices = find(Partition.grid_x(:,3) == Height);

if isempty(Indices)
    error('The required initial condition is not compatible with the state partition');
else
    gridXY = Partition.grid_x(:,1:2);
    x = unique(gridXY(:,1));
    y = unique(gridXY(:,2));
    [X,Y] = meshgrid(x,y);
    
    tempValue = ValueFunction(Indices);
    tempPolicy = OptInput(Indices);
    
    ProjValueFunction = reshape(tempValue,size(X));
    OptPolicy = reshape(tempPolicy,size(X));
end


end

% 
% [XX,YY,ProjValueFunc,OptPolicy] = plot_results(InitialState_X3,valueFunction(:,1),OptInput(:,1),Grid.getValues,param);
% 
% levels = [0.1,0.2,0.4,0.5,0.6,0.8,0.85,0.9,1];
% 
% figure 
% h = contour(XX,YY,ProjValueFunc,levels,'ShowText','on','LineWidth',2);
% 
% 
% figure
% pcolor(XX,YY,OptPolicy)


% 
% InitialState_X3 = 420;
% InitialState_X3 = Grid.X3(1,1,find(Grid.X3(1,1,:) <= InitialState_X3,1,'last'));
% 
% [XX,YY,ProjValueFunc,OptPolicy] = plot_results(InitialState_X3,valueFunction(:,1),OptInput(:,1),Grid,param);
% 
% levels = [0.1,0.2,0.4,0.5,0.6,0.8,0.85,0.9,1];
% 
% figure
% h = contour(XX,YY,ProjValueFunc,levels,'ShowText','on','LineWidth',2);
% 
% 
% figure
% pcolor(XX,YY,OptPolicy)




% 
% InitialState_X3 = 600;
% InitialState_X3 = Grid.X3(1,1,find(Grid.X3(1,1,:) <= InitialState_X3,1,'last'));
% 
% [XX,YY,ProjValueFunc,OptPolicy] = plot_results(InitialState_X3,valueFunction(:,1),OptInput(:,1),Grid,param);
% 
% figure 
% h = contour(XX,YY,ProjValueFunc,levels,'ShowText','on','LineWidth',2);

%Functions related to the partition and generating an estimate of the transition kernel




