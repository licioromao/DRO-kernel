function out = getDateSaveFile(TimeHorizon,NumberOfPartitions,NumberOfMonteCarlo,TypeOfVectorField,param)

ep = []; rhoMu = []; rhoSigma = [];

if isfield(param,'ep')
    ep = param.ep;
end

if isfield(param,'rhoMu')
    rhoMu = param.rhoMu;
end

if isfield(param,'rhoSigma')
    rhoSigma = param.rhoSigma;
end

DestinationFolder = sprintf('./Results/results_%s',char(java.net.InetAddress.getLocalHost.getHostName));
[status{1},msg{1}] = mkdir(DestinationFolder);

switch TypeOfVectorField
    case 'TCL'
        DestinationFolder1 = strcat(sprintf('%s/',DestinationFolder),TypeOfVectorField);
        [status{2},msg{2}] = mkdir(DestinationFolder1);
        
        DestinationFolder2 = strcat('/N_',sprintf('%d',TimeHorizon),'_Partition_',sprintf('%d_',NumberOfPartitions),'m_',sprintf('%d',NumberOfMonteCarlo));
        DirectoryPath = strcat(DestinationFolder1,DestinationFolder2);
        [status{3},msg{3}] = mkdir(DirectoryPath);
     case 'ChainInt'
        DestinationFolder1 = strcat(sprintf('%s/',DestinationFolder),TypeOfVectorField);
        [status{2},msg{2}] = mkdir(DestinationFolder1);
        
        DestinationFolder2 = strcat('/N_',sprintf('%d',TimeHorizon),'_Partition_',sprintf('%d_%d_',NumberOfPartitions),'m_',sprintf('%d',NumberOfMonteCarlo));
        DirectoryPath = strcat(DestinationFolder1,DestinationFolder2);
        [status{3},msg{3}] = mkdir(DirectoryPath);
        
    case 'Fishery'
        DestinationFolder1 = strcat(sprintf('%s/',DestinationFolder),TypeOfVectorField);
        [status{2},msg{2}] = mkdir(DestinationFolder1);
        
        DestinationFolder2 = strcat('/N_',sprintf('%d',TimeHorizon),'_Partition_',sprintf('%d_%d_%d_',NumberOfPartitions),'m_',sprintf('%d',NumberOfMonteCarlo));
        DirectoryPath = strcat(DestinationFolder1,DestinationFolder2);
        [status{3},msg{3}] = mkdir(DirectoryPath);
        
    case 'CarPole'
        DestinationFolder1 = strcat(sprintf('%s/',DestinationFolder),TypeOfVectorField);
        [status{2},msg{2}] = mkdir(DestinationFolder1);
        
        DestinationFolder2 = strcat('/N_',sprintf('%d',TimeHorizon),'_Partition_',sprintf('%d_%d_%d_%d_',NumberOfPartitions),'m_',sprintf('%d',NumberOfMonteCarlo));
        DirectoryPath = strcat(DestinationFolder1,DestinationFolder2);
        [status{3},msg{3}] = mkdir(DirectoryPath);
        
    case 'CarPoleNL'
        DestinationFolder1 = strcat(sprintf('%s/',DestinationFolder),TypeOfVectorField);
        [status{2},msg{2}] = mkdir(DestinationFolder1);
        
        DestinationFolder2 = strcat('/N_',sprintf('%d',TimeHorizon),'_Partition_',sprintf('%d_%d_%d_%d_',NumberOfPartitions),'m_',sprintf('%d',NumberOfMonteCarlo));
        DirectoryPath = strcat(DestinationFolder1,DestinationFolder2);
        [status{3},msg{3}] = mkdir(DirectoryPath);
        
    otherwise
        NotImplemented();
end


if ~isempty(ep)
    AppendPath = sprintf('/ep_%.4f',ep);
    DirectoryPath = strcat(DirectoryPath,AppendPath);
    [status{4},msg{4}] = mkdir(DirectoryPath);
    if ~isempty(rhoMu)
        AppendPath = sprintf('/rhoMu_%.2f',rhoMu);
        DirectoryPath = strcat(DirectoryPath,AppendPath);
        [status{5},msg{5}] = mkdir(DirectoryPath);
        if ~isempty(rhoSigma)
            AppendPath = sprintf('/rhoSigma_%.2f',rhoSigma);
            DirectoryPath = strcat(DirectoryPath,AppendPath);
            [status{6},msg{6}] = mkdir(DirectoryPath);
            
            if ~(status{1}&& status{2} &&status{3} && status{4} &&status{5} && status{6} )
                error('%s\n%s\n%s\n%s\n%s\n%s\n',msg{1},msg{2},msg{3},msg{4},msg{5},msg{6});
            end
        end
    end
elseif ~isempty(rhoMu)
    AppendPath = sprintf('/rhoMu_%.2f',rhoMu);
    DirectoryPath = strcat(DirectoryPath,AppendPath);
    [status{4},msg{4}] = mkdir(DirectoryPath);
    if ~isempty(rhoSigma)
        AppendPath = sprintf('/rhoSigma_%.2f',rhoSigma);
        DirectoryPath = strcat(DirectoryPath,AppendPath);
        
        [status{5},msg{5}] = mkdir(DirectoryPath);
        
        if ~(status{1}&& status{2} &&status{3} && status{4} &&status{5} )
            error('%s\n%s\n%s\n%s\n%s\n',msg{1},msg{2},msg{3},msg{4},msg{5});
        end
        
    end
end


DayTime = clock;
Months{1} = 'Jan'; Months{2} = 'Feb'; Months{3} = 'Mar'; Months{4} = 'Apr';
Months{5} = 'May'; Months{6} = 'Jun'; Months{7} = 'Jul'; Months{8} = 'Aug';
Months{9} = 'Sep'; Months{10} = 'Oct'; Months{11} = 'Nov'; Months{12} = 'Dec';

YY = sprintf('%d',DayTime(1)); MM = Months{DayTime(2)}; DD = sprintf('%d',DayTime(3));
HH = sprintf('%d',DayTime(4)); Min = sprintf('%d',DayTime(5)); Sec = sprintf('%d',floor(DayTime(6)));

FileName = strcat('Date_',DD,'_',MM,'_',YY,'_',HH,'_',Min,'_',Sec);

out.FullPath = strcat(DirectoryPath,'/',FileName,'.mat');

out.Directory = DirectoryPath;
out.FileName = FileName;

end

