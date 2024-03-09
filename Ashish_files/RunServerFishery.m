function out = RunServerFishery(N,NumberOfPartitions,NumberMonteCarlo)

AmbiguityTypes = {'NoAmbiguity','KernelAmbiguity'};

out = Fishery_management(int16(N),NumberOfPartitions,NumberMonteCarlo,AmbiguityTypes);

end



