%%% https://www.briandalessandro.com/blog/how-to-make-a-borderless-subplot-of-images-in-matlab/
function h = subplottight(n,m,i)
    [c,r] = ind2sub([m n], i);
    ax = subplot('Position', [(c-1)/m, 1-(r)/n, 1/m, 1/n])
    set(ax,'xtick',[])
    set(ax,'xticklabel',[])
    set(ax,'ytick',[])
    set(ax,'yticklabel',[])
    if(nargout > 0)
      h = ax;
    end
end