function val=catch_22_feeder(data)
        [rows,~]=size(data);
        vals=[];
        for i =1:rows
            label=data(i,1);
            x=data(i,2:end);
            lastidx = find(isnan(x),1) - 1;
            x=x(1:lastidx)';
            [new,~]=catch22_all(x);
            l=[label,new'];
            vals=[vals;l];    
        end
        val=vals;
end