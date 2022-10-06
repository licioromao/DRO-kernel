function StringXU = createXandUString(x,u)

tempX = compose('%.4f',x);

if ~isempty(u)
    tempU = compose('%.4f',u);
    StringXU = strcat('X:',tempX{:},'|','U:',tempU{:});
else
    StringXU = strcat('X:',tempX{:});
end

end

