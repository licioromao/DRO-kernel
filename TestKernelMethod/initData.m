% Adding classes definitions that are necessary to import the data in
% TestKernel.mat

addpath '/Users/licioromao/Dropbox/Collaborations/Ashish and Licio/Distributionally robust control problems/Licio_code/AmbiguityClasses' 
addpath '/Users/licioromao/Dropbox/Collaborations/Ashish and Licio/Distributionally robust control problems/Licio_code/VectorFields'
addpath '/Users/licioromao/Dropbox/Collaborations/Ashish and Licio/Distributionally robust control problems/Licio_code/PartitionClassesAndFunctions'
addpath '/Users/licioromao/Dropbox/Collaborations/Ashish and Licio/Distributionally robust control problems/Licio_code/ValueFunctionComputation'
addpath '/Users/licioromao/Dropbox/Collaborations/Ashish and Licio/Distributionally robust control problems/Licio_code/AuxiliaryFunctions'
addpath '/Users/licioromao/Dropbox/Collaborations/Ashish and Licio/Distributionally robust control problems/Licio_code/Kernels'


load('TestKernel.mat','OptValueFunc','OptPolicy','Grid','param')

gridX = Grid.getValues.Partition.grid_x;