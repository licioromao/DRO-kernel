function out = get_date_save_file(time_horizon,number_of_points,number_of_MC_simulations,...
                                type_vector_field,param)

ep = []; radius_mean = []; radius_variance = [];

if isfield(param,'ep')
    ep = param.ep;
end

if isfield(param,'radius_mean')
    radius_mean = param.radius_mean;
end

if isfield(param,'radius_variance')
    radius_variance = param.radius_variance;
end


destination_folder = strcat(param.path_project,sprintf('/Results/results_%s',char(java.net.InetAddress.getLocalHost.getHostName)));
[status{1},msg{1}] = mkdir(destination_folder);

switch type_vector_field
    case 'TCL'
        destination_folder_1 = strcat(sprintf('%s/',destination_folder),type_vector_field);
        [status{2},msg{2}] = mkdir(destination_folder_1);
        
        destination_folder_2 = strcat('/N_',sprintf('%d',time_horizon),'_Partition_',sprintf('%d_',number_of_points),'m_',sprintf('%d',number_of_MC_simulations));
        directory_path = strcat(destination_folder_1,destination_folder_2);
        [status{3},msg{3}] = mkdir(directory_path);
     case 'ChainInt'
        destination_folder_1 = strcat(sprintf('%s/',destination_folder),type_vector_field);
        [status{2},msg{2}] = mkdir(destination_folder_1);
        
        destination_folder_2 = strcat('/N_',sprintf('%d',time_horizon),'_Partition_',sprintf('%d_%d_',number_of_points),'m_',sprintf('%d',number_of_MC_simulations));
        directory_path = strcat(destination_folder_1,destination_folder_2);
        [status{3},msg{3}] = mkdir(directory_path);
        
    case 'Fishery'
        destination_folder_1 = strcat(sprintf('%s/',destination_folder),type_vector_field);
        [status{2},msg{2}] = mkdir(destination_folder_1);
        
        destination_folder_2 = strcat('/N_',sprintf('%d',time_horizon),'_Partition_',sprintf('%d_%d_%d_',number_of_points),'m_',sprintf('%d',number_of_MC_simulations));
        directory_path = strcat(destination_folder_1,destination_folder_2);
        [status{3},msg{3}] = mkdir(directory_path);
        
    case 'CarPole'
        destination_folder_1 = strcat(sprintf('%s/',destination_folder),type_vector_field);
        [status{2},msg{2}] = mkdir(destination_folder_1);
        
        destination_folder_2 = strcat('/N_',sprintf('%d',time_horizon),'_Partition_',sprintf('%d_%d_%d_%d_',number_of_points),'m_',sprintf('%d',number_of_MC_simulations));
        directory_path = strcat(destination_folder_1,destination_folder_2);
        [status{3},msg{3}] = mkdir(directory_path);
        
    case 'CarPoleNL'
        destination_folder_1 = strcat(sprintf('%s/',destination_folder),type_vector_field);
        [status{2},msg{2}] = mkdir(destination_folder_1);
        
        destination_folder_2 = strcat('/N_',sprintf('%d',time_horizon),'_Partition_',sprintf('%d_%d_%d_%d_',number_of_points),'m_',sprintf('%d',number_of_MC_simulations));
        directory_path = strcat(destination_folder_1,destination_folder_2);
        [status{3},msg{3}] = mkdir(directory_path);
        
    otherwise
        not_implemented();
end


if ~isempty(ep)
    append_path = sprintf('/ep_%.4f',ep);
    directory_path = strcat(directory_path,append_path);
    [status{4},msg{4}] = mkdir(directory_path);
    if ~isempty(radius_mean)
        append_path = sprintf('/rhoMu_%.2f',radius_mean);
        directory_path = strcat(directory_path,append_path);
        [status{5},msg{5}] = mkdir(directory_path);
        if ~isempty(radius_variance)
            append_path = sprintf('/rhoSigma_%.2f',radius_variance);
            directory_path = strcat(directory_path,append_path);
            [status{6},msg{6}] = mkdir(directory_path);
            
            if ~(status{1}&& status{2} &&status{3} && status{4} &&status{5} && status{6} )
                error('%s\n%s\n%s\n%s\n%s\n%s\n',msg{1},msg{2},msg{3},msg{4},msg{5},msg{6});
            end
        end
    end
elseif ~isempty(radius_mean)
    append_path = sprintf('/rhoMu_%.2f',radius_mean);
    directory_path = strcat(directory_path,append_path);
    [status{4},msg{4}] = mkdir(directory_path);
    if ~isempty(radius_variance)
        append_path = sprintf('/rhoSigma_%.2f',radius_variance);
        directory_path = strcat(directory_path,append_path);
        
        [status{5},msg{5}] = mkdir(directory_path);
        
        if ~(status{1}&& status{2} &&status{3} && status{4} &&status{5} )
            error('%s\n%s\n%s\n%s\n%s\n',msg{1},msg{2},msg{3},msg{4},msg{5});
        end
        
    end
end


day_time = clock;
months{1} = 'Jan'; months{2} = 'Feb'; months{3} = 'Mar'; months{4} = 'Apr';
months{5} = 'May'; months{6} = 'Jun'; months{7} = 'Jul'; months{8} = 'Aug';
months{9} = 'Sep'; months{10} = 'Oct'; months{11} = 'Nov'; months{12} = 'Dec';

YY = sprintf('%d',day_time(1)); MM = months{day_time(2)}; DD = sprintf('%d',day_time(3));
HH = sprintf('%d',day_time(4)); Min = sprintf('%d',day_time(5)); Sec = sprintf('%d',floor(day_time(6)));

file_name = strcat('Date_',DD,'_',MM,'_',YY,'_',HH,'_',Min,'_',Sec);

out.full_path = strcat(directory_path,'/',file_name,'.mat');

out.directory = directory_path;
out.file_name = file_name;

end

