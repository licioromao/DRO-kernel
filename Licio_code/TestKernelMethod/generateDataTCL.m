function [DATA,PartitionState] = generateDataTCL(Grid,OptPolicy,m,param)


gridX = Grid.getValues.Partition.grid_x;

x = linspace(gridX(1),gridX(end),m)';
u = zeros(length(x),1);
y = zeros(length(x),1);

for i =1:length(x)
    index = Grid.getElementPartition(x(i)).index;
    u(i) = OptPolicy(index);
    tempVF = VectorFieldTCL(x(i),u(i),param);
    tempVF.Noise = generateNoise(param,'TCL');
    y(i) = tempVF.IterateDynamics();
end

DATA = [y,x,u];
PartitionState = StatePartition(m,param.SafeSet,'TCL');

end