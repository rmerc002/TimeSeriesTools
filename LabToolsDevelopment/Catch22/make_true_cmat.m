function val=make_true_cmat(x,y)
mat_lim=unique(y);
lims=size(mat_lim);
conf_mat=zeros(lims(1),lims(1));
loop_s=size(x);
for i=1:loop_s(1)
    if x(i)==y(i)
        a=find(mat_lim==x(i));
        conf_mat(a,a)=conf_mat(a,a)+1;
    else 
        a=find(mat_lim==x(i));
        b=find(mat_lim==y(i));
        conf_mat(a,b)=conf_mat(a,b)+1;
    end
end
sumu=sum(conf_mat,2);
%conf_mat=conf_mat./sumu;
val=conf_mat;
end