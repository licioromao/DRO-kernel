function out = getDateSaveFile(NumberOfPartitions,TypeOfVectorField)

switch TypeOfVectorField
    case 'Fishery'
        FirstPart = sprintf('%s_Partition_%d_%d_%d_Date_',TypeOfVectorField,NumberOfPartitions);
    case 'TCL'
        FirstPart = sprintf('%s_Partition_%d_Date_',TypeOfVectorField,NumberOfPartitions);
    otherwise
end


DayTime = clock;
Months{1} = 'Jan'; Months{2} = 'Feb'; Months{3} = 'Mar'; Months{4} = 'Apr';
Months{5} = 'May'; Months{6} = 'Jun'; Months{7} = 'Jul'; Months{8} = 'Aug';
Months{9} = 'Sep'; Months{10} = 'Oct'; Months{11} = 'Nov'; Months{12} = 'Dec';

YY = sprintf('%d',DayTime(1)); MM = Months{DayTime(2)}; DD = sprintf('%d',DayTime(3));
HH = sprintf('%d',DayTime(4)); Min = sprintf('%d',DayTime(5)); Sec = sprintf('%d',floor(DayTime(6)));
out = strcat('./results/',FirstPart,DD,'_',MM,'_',YY,'_',HH,'_',Min,'_',Sec);

end

