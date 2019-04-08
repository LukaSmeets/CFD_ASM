function grid = meshgrid(x, NPI)

grid = [];
    for i = 1:length(x)/NPI 
        grid = [grid, x((i-1)*NPI+1:NPI*i)];
    end
    
end