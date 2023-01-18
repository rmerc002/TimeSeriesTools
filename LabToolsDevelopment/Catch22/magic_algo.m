function [val,prob]=magic_algo(x,y,labels_train)
prob=[];
k=1;


    v=pdist2(y,x,'euclidean');
    if k==1

    [r,index]=min(v,[],2); 
      
    
    else
    [r,index]=mink(v,k,1);
    end


val=index;
prob=r;
v=[];
end