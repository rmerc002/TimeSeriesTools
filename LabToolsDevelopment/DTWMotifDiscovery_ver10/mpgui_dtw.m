classdef mpgui_dtw < handle
    properties
        fig;
        data;
        sublen;
        datalen;
        proflen;
        profylims;
        dataylims;
        dataAx;
        profileAx;
        motifAxed;
        motifAxdtw;        
        vis;
        stopBtn;
        shouldHalt;
    end
    
    properties(Constant, Access = private)
        discordColor = 'brg';
        motifColor = 'gc';
        neighborColor = [0.5, 0.5, 0.5];
    end
    
    methods
        
        function gui = mpgui_dtw(data, sublen)            
            gui.data = data;
            gui.sublen = sublen;
            gui.datalen = length(data);
            gui.proflen = gui.datalen - sublen + 1;
            gui.profylims = [0, floor(2 * sqrt(sublen))];
            gui.dataylims = ([min(data), max(data)]);
            sh = @(fig, obj) gui.stopCallback();
            gui.fig = figure('visible', 'off', 'toolbar', 'none', 'Name', 'SWAMP');
            gui.dataAx = subplot('Position', [0.1 0.8 0.75 0.15]);
            gui.profileAx = subplot('Position', [0.1 0.55 0.75 0.15]);
            gui.motifAxed = subplot('Position', [0.1 0.3 0.75 0.15]);
            gui.motifAxdtw = subplot('Position', [0.1 0.05 0.75 0.15]);
            gui.stopBtn = uicontrol('parent', gui.fig, 'style', 'pushbutton', 'units', 'normalized', 'Position', [0.88 0.05 0.1 0.075], 'string', 'Stop', 'callback', sh);
            gui.shouldHalt = false;
            gui.vis = false;
        end
        
        function delete(gui)
            % very basic destructor, it just ensures that gui.close is
            % called when the object is deleted if it was never revealed.
            % It's better to close these manually rather than rely on this
            if ~gui.vis
                gui.close();
            end
        end
        
        function drawgui(gui)
            if ~gui.vis
                gui.vis = true;
                gui.fig.Visible = 'on';
            end
            drawnow;
        end
        
        function plotProfile(gui, prof)
            plot(gui.profileAx, 1 : gui.proflen, prof, 'k');
            xlim(gui.profileAx, [1 gui.datalen]);
            xticks(gui.profileAx, [1, gui.datalen]);
            ylim(gui.profileAx, gui.profylims);
            yticks(gui.profileAx, gui.profylims);
            ytickformat(gui.profileAx, '%.4f');
            t = title(gui.profileAx, 'The Euclidean distance Matrix Profile', 'horizontalalignment', 'left', 'units', 'normalized', 'fontweight', 'normal');
            p = get(t, 'position');
            set(t, 'position', [0 p(2) p(3)]);
        end
        
        function plotData(gui, motifidx, neighboridx)
            plot(gui.dataAx, 1 : gui.datalen, gui.data, 'Color' , [0.75 0.75 0.75]);
            if nargin == 3
                hold(gui.dataAx, 'on');
                plot(gui.dataAx, motifidx : motifidx + gui.sublen - 1, gui.data(motifidx : motifidx + gui.sublen-1), gui.motifColor(1), 'Linewidth', 2);
                plot(gui.dataAx, neighboridx : neighboridx + gui.sublen - 1, gui.data(neighboridx : neighboridx + gui.sublen - 1), gui.motifColor(2), 'Linewidth', 2);
                hold(gui.dataAx, 'off');
                
                t = title(gui.dataAx, 'The input time series. Motifs are color coded.', 'horizontalalignment', 'left', 'units', 'normalized', 'fontweight', 'normal');
                %t = title(gui.dataAx, sprintf('%.1f%% done: The input time series. Motifs are color coded', progress),'horizontalalignment', 'left', 'units', 'normalized', 'fontweight', 'normal');
            else
                t = title(gui.dataAx, 'The input time series', 'horizontalalignment', 'left', 'units', 'normalized', 'fontweight', 'normal');
            end
            xlim(gui.dataAx, [1 gui.datalen]);
            xticks(gui.dataAx, [1 gui.datalen]);
            ylim(gui.dataAx, gui.dataylims);
            yticks(gui.dataAx, gui.dataylims);
            ytickformat(gui.dataAx, '%.2f');
            p = get(t, 'position');
            set(t, 'position', [0 p(2) p(3)]);
        end
        
        function plotMotif(gui, axhandle, meas, motifidx, neighboridx)
            
            mot = (gui.data(motifidx : motifidx + gui.sublen - 1)...
                - mean(gui.data(motifidx : motifidx + gui.sublen - 1)))...
                ./std(gui.data(motifidx : motifidx + gui.sublen - 1), 1);
            neigh = (gui.data(neighboridx : neighboridx + gui.sublen - 1)...
                   - mean(gui.data(neighboridx : neighboridx + gui.sublen - 1)))...
                   ./std(gui.data(neighboridx : neighboridx + gui.sublen - 1), 1);
            plot(axhandle, 1 : gui.sublen, mot, gui.motifColor(1));
            hold(axhandle, 'on');
            plot(axhandle, 1 : gui.sublen, neigh, gui.motifColor(2));
            xlim(axhandle, [1 gui.sublen]);
            y = [floor(min(mot)) ceil(max(mot))];
            ylim(axhandle, y);
            yticks(axhandle, y);
            ytickformat(axhandle, '%.2f');
            t = title(axhandle, sprintf('The top 1 motif pair using %s is located at %d and %d', meas, motifidx, neighboridx), 'horizontalalignment', 'left', 'units', 'normalized', 'fontweight', 'normal');
            p = get(t, 'position');
            set(t, 'position', [0 p(2) p(3)]);
            hold(axhandle, 'off');
        end
        
        function plotMotifed(gui, motifidx, neighboridx)
            gui.plotMotif(gui.motifAxed, 'euclidean distance', motifidx, neighboridx);
            
        end
        
        function plotMotifdtw(gui, motifidx, neighboridx)
            gui.plotMotif(gui.motifAxdtw, 'dynamic time warping', motifidx, neighboridx);
            %hold(axhandle, 'off');
        end
        
        function stopCallback(gui)
            gui.shouldHalt = true;                        
            drawnow;
        end
        
        function close(gui)
            if isvalid(gui.fig)
                close(gui.fig);
            end
        end
        
    end
end
