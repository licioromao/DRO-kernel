function L = PlotResults(FILE)

load(FILE,'-regexp','ValueFunc*');
load(FILE,'Grid');

L = size(whos,1) - 2;

switch FILE(11:13)
    case 'TCL'
        X = Grid.getValues.Partition.X
        tempVariable(:) = whos('-regexp','ValueFunc*');   
        figure;
        hold on
        Name = cell(L,1);
        for i =1:L
            grid on;
            box on;
            Name{i} = tempVariable(1,i).name;
            plot(X,eval(Name{i}).ValueFunction(1:end-1,1),'LineWidth',1.5);
        end
        legend(Name)
    case 'Fis'
    otherwise
        NotImplemented();
end

end

