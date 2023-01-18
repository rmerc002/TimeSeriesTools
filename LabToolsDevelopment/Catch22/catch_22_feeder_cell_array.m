function output=catch_22_feeder_cell_array(data)
        [rows,~] = size(data);
        output = [];
        for i =1:rows
            [new,~] = catch22_all(data{i,1});
            output = [output;new'];    
        end
end