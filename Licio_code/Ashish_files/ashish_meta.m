%run('testTCL.m');

FisheryPartition = [3,3,3];
mFishery = 5;
ep = 0.01;
rhoMu = [10,20,50];
rhoSigma = [5,10,20];
run('testFishery.m');
