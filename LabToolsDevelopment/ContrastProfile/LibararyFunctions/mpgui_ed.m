classdef mpgui_ed < handle
    properties
        fig;
        data;
        subseqlen;
        datalen;
        subcount;
        profylims;
        dataylims;
        dataAx;
        profileAx;
        motifAx1;
        motifAx2;
        motifAx3;
        discordAx;
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
        
        function gui = mpgui_ed(data, subseqlen, is_interactive)
            if nargin == 2
                is_interactive = false;
            end
            gui.data = data;
            gui.subseqlen = subseqlen;
            gui.datalen = length(data);
            gui.subcount = gui.datalen - subseqlen + 1;
            gui.profylims = [0, 2 * sqrt(subseqlen)];
            gui.dataylims = 1.05 * [min(data), max(data)];
            sh = @(fig, obj) gui.stopCallback();
            gui.fig = figure('visible', 'off', 'toolbar', 'none');
            height = 0.08;
            top = 0.85;
            left = 0.08;
            vdiff = 0.16;
            width = 0.75;
            gui.dataAx = subplot('Position', [left top width height]);
            gui.profileAx = subplot('Position',[left top-vdiff   width height]);
            gui.motifAx1 = subplot('Position', [left top-2*vdiff width height]);
            gui.motifAx2 = subplot('Position', [left top-3*vdiff width height]);
            gui.motifAx3 = subplot('Position', [left top-4*vdiff width height]);
            gui.discordAx = subplot('Position',[left top-5*vdiff width height]);
            if is_interactive
                gui.stopBtn = uicontrol('parent', gui.fig, 'style', 'pushbutton', 'units', 'normalized', 'Position', [0.88 0.05 0.1 0.075], 'string', 'Stop', 'callback', sh);
                gui.shouldHalt = false;
            end
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
            plot(gui.profileAx, 1 : gui.subcount, prof, 'b');
            xlim(gui.profileAx, [1 gui.subcount]);
            xticks(gui.profileAx, [1, gui.subcount]);
            ylim(gui.profileAx, gui.profylims);
            yticks(gui.profileAx, gui.profylims);
            ytickformat(gui.profileAx, '%.0f');
            t = title(gui.profileAx, 'The best-so-far matrix profile', 'horizontalalignment', 'left', 'units', 'normalized', 'fontweight', 'normal');
            p = get(t, 'position');
            set(t, 'position', [0 p(2) p(3)]);
        end
        
        function plotData(gui, motifidx, neighboridx)
            % second and third arg are optional
            plot(gui.dataAx, 1 : gui.datalen, gui.data, 'b');
            if nargin == 3
                hold(gui.dataAx, 'on');
                plot(gui.dataAx, motifidx : motifidx + gui.subseqlen - 1, gui.data(motifidx : motifidx + gui.subseqlen-1), gui.motifColor(1));
                plot(gui.dataAx, neighboridx : neighboridx + gui.subseqlen - 1, gui.data(neighboridx : neighboridx + gui.subseqlen - 1), gui.motifColor(2));
                hold(gui.dataAx, 'off');
                t = title(gui.dataAx, 'The input time series, motifs are color coded', 'horizontalalignment', 'left', 'units', 'normalized', 'fontweight', 'normal');
            else
                t = title(gui.dataAx, 'The input time series', 'horizontalalignment', 'left', 'units', 'normalized', 'fontweight', 'normal');
            end
            xlim(gui.dataAx, [1 gui.datalen]);
            xticks(gui.dataAx, [1 gui.datalen]);
            ylim(gui.dataAx, gui.dataylims);
            yticks(gui.dataAx, []);
            set(gui.dataAx, 'YColor', [1 1 1])
            p = get(t, 'position');
            set(t, 'position', [0 p(2) p(3)]);
        end
        
        function plotMotif(gui, axhandle, motifidx, neighboridx, motiftitle)
            
            mot = (gui.data(motifidx : motifidx + gui.subseqlen - 1)...
                - mean(gui.data(motifidx : motifidx + gui.subseqlen - 1)))...
                ./std(gui.data(motifidx : motifidx + gui.subseqlen - 1), 1);
            % mild hack 1/14/21 since there are too many cases to deal with
            % when missing data is encountered
            % The approaches in the very early (predating mpx) codes
            % are blatantly incorrect in some cases
            % and should not be used as reference.
            if any(isnan(mot))
                warning('non-normalizable motif sequence');
                return
            end
            neigh = (gui.data(neighboridx : neighboridx + gui.subseqlen - 1)...
                   - mean(gui.data(neighboridx : neighboridx + gui.subseqlen - 1)))...
                   ./std(gui.data(neighboridx : neighboridx + gui.subseqlen - 1), 1);
            plot(axhandle, 1 : gui.subseqlen, mot, gui.motifColor(1));
            hold(axhandle, 'on');
            plot(axhandle, 1 : gui.subseqlen, neigh, gui.motifColor(2));
            xlim(axhandle, [1 gui.subseqlen]);
            y = [min(mot) max(mot)];
            ylim(axhandle, y);
            yticks(axhandle, []);
            set(axhandle, 'YColor', [1 1 1])
            t = title(axhandle, motiftitle, 'horizontalalignment', 'left', 'units', 'normalized', 'fontweight', 'normal');
            p = get(t, 'position');
            set(t, 'position', [0 p(2) p(3)]);
        end
        
        function plotDiscords(gui, discordIdxs, discordtext)
            disclims = zeros(2, 1);
            disc = gui.data(discordIdxs(1) : discordIdxs(1) + gui.subseqlen - 1);
            disc = (disc - mean(disc)) / std(disc, 1);
            if any(isnan(disc))
                % mild hack to deal with large amounts of missing data
                % without a crash. Don't refer to motif and discord 
                % handling embedded in methods
                % predating mpx here. Those are completely inconsisent
                % and wrong.
                % edit: 1/14/21
                warning('non-normalizable discord encountered');
                return
            end
            disclims = [min(disclims(1), min(disc)); max(disclims(2), max(disc))];
            plot(gui.discordAx, disc, gui.discordColor(1));
            hold(gui.discordAx, 'on');
            for j = 2 : length(discordIdxs)
                disc = gui.data(discordIdxs(j) : discordIdxs(j) + gui.subseqlen - 1);
                disc = (disc - mean(disc)) / std(disc, 1);
                plot(gui.discordAx, disc, gui.discordColor(j));
                disclims = [min(disclims(1), min(disc)); max(disclims(2), max(disc))];
            end
            hold(gui.discordAx, 'off');
            t = title(gui.discordAx, discordtext, 'horizontalalignment', 'left', 'units', 'normalized', 'fontweight', 'normal');
            xlim(gui.discordAx, [1 gui.subseqlen]);
            ylim(gui.discordAx, disclims);
            yticks(gui.discordAx, []);
            set(gui.discordAx, 'YColor', [1 1 1])
            p = get(t, 'position');
            set(t, 'position', [0 p(2) p(3)]);
            if disclims(1) ~= disclims(2)
                ylim(gui.discordAx, 1.05 * disclims);
            end
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
