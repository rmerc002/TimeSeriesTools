function RI = RandIndex(cls1, cls2)

    L=length(cls1);
    A=0;B=0;C=0;D=0;
    for i=1:L
        for j=i+1:L
            
            if((cls1(i)==cls1(j))&&(cls2(i)==cls2(j)))
                A=A+1;
            elseif((cls1(i)~=cls1(j))&&(cls2(i)~=cls2(j)))
                B=B+1;
            elseif((cls1(i)==cls1(j))&&(cls2(i)~=cls2(j)))
                C=C+1;
            else
                D=D+1;
            end
        end
    end
    
    RI = (A+B)/(A+B+C+D);

end