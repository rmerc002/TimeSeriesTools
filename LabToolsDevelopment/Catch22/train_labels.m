function [val,prob]=train_labels(x)
prob=[];
k=1;
u_size=size(x);

    v=pdist2(x,x,'euclidean');
    
    v_iden=eye(u_size(1));
    v_iden=reshape(v_iden,[],1);
    v_iden(v_iden==1)=99999;
    v_iden=reshape(v_iden,[u_size(1),u_size(1)]);
    v_final=v+v_iden;
    if k==1

    [r,index]=min(v_final,[],2); 
    else
    [r,index]=mink(v_final,k,1);
    end
v=[];
val=index;
end