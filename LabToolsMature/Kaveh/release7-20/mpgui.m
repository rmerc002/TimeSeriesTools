classdef mpgui < handle
    
    properties(Access = private)
        fig;
        data;
        mu;
        invnorm;
        subLen;
        dataAx;
        profileAx;
        discordAx;
        motifAx;
        motifNeighbLineH;
        dataText;
        profileText;
        discordText;
        motifText;
        discardBtn;
        discards;
        stopBtn;
        shouldHalt_;
    end
    
    properties(Constant, Access = private)
        txtTemp = {
            'The best motif pair is located at %d (green) and %d (cyan)';
            'The 2nd best motif pair is located at %d (green) and %d (cyan)';
            'The 3rd best motif pair is located at %d (green) and %d (cyan)';
            };
        discordColor = 'brg';
        motifColor = 'gc';
        neighborColor = [0.5, 0.5, 0.5];
    end
    
    methods
        function setTitle(gui, txt)
            gui.fig.Name = txt;
            drawnow;
        end
    end
    
    methods(Access = private)
        
        function gui = mpgui(data, mu, invnorm, subLen)
            dataLen = length(data);
            gui.data = data;
            gui.mu = mu;
            gui.invnorm = invnorm;
            gui.subLen = subLen;
            gui.discardBtn = gobjects(3, 1);
            rs = @(fig, obj) gui.mainResize();
            sh = @(fig, obj) gui.stopCallback();
            gui.fig = figure('visible', 'off', 'toolbar', 'none', 'ResizeFcn', rs);
            lims = [-0.05 1.05];
            gui.dataAx = axes(gui.fig, 'units', 'pixels', 'xlim', [1, dataLen], 'xtick', [1, dataLen], 'ylim', 1.05 * [min(data); max(data)], 'ytick', [], 'ycolor', [1, 1, 1]);
            gui.profileAx = axes(gui.fig, 'units', 'pixels', 'xlim', [1, dataLen], 'xtick', [1, dataLen], 'ylim', [0, 2 * sqrt(subLen)]);
            gui.discordAx = axes(gui.fig, 'units', 'pixels', 'xlim', [1, subLen], 'xtick', [1, subLen], 'ylim', lims, 'ytick', [], 'ycolor', [1, 1, 1]);
            gui.stopBtn = uicontrol('parent', gui.fig, 'style', 'pushbutton', 'string', 'Stop', 'fontsize', 10, 'callback', sh);
            gui.dataText = uicontrol('parent', gui.fig, 'style', 'text', 'string', '', 'fontsize', 10, 'horizontalalignment', 'left');
            gui.profileText = uicontrol('parent', gui.fig, 'style', 'text', 'string', 'The best-so-far matrix profile', 'horizontalalignment', 'left', 'fontsize', 10);
            gui.discordText = uicontrol('parent', gui.fig, 'style', 'text', 'string', '', 'horizontalalignment', 'left', 'fontsize', 10);
            gui.motifAx = gobjects(3, 1);
            for i = 1:3
                gui.motifAx(i) = axes('parent', gui.fig, 'units', 'pixels', 'xlim', [1, subLen], 'xtick', [1, subLen], 'ylim', lims, 'ytick', [], 'ycolor', [1, 1, 1]);
            end
            gui.motifText = gobjects(3,1);
            for i = 1:3
                gui.motifText(i) = uicontrol('parent', gui.fig, 'style', 'text', 'string', '', 'horizontalalignment', 'left', 'fontsize', 10);
            end
            figPosition = get(gui.fig, 'Position');
            set(gui.stopBtn, 'position', [figPosition(3) - 120, 30, 90, 20]);
            gui.shouldHalt_ = false;
        end
        
        % The very old gui normalized explicitly to the 0 to 1 range. That was
        % didn't really reflect a shape based comparison. This one relies
        % on whatever normalization parameters were passed in when the gui
        % is initialized.
        
        function ss = normalizedSubseq(gui, idx)
            ss = (gui.data(idx : idx + gui.subLen - 1) - gui.mu(idx)) .* gui.invnorm(idx);
        end
        
        function discard(gui)
            gui.discardBtn = gobjects(3);
            for i = 1:3
                gui.discardBtn(i) = uicontrol('parent', gui.fig, 'style', 'pushbutton', 'string', 'Discard', 'fontsize', 10, 'callback', @(src, cbdata) pushDiscardBtn(src, cbdata, i));
            end
        end
        
        function plotData(gui)
            hold(gui.dataAx, 'on');
            plot(gui.dataAx, gui.data, 'r');
            hold(gui.dataAx, 'off');
        end
        
        function plotProfile(gui, matrixProfile)
            hold(gui.profileAx, 'on');
            plot(gui.profileAx, 1 : length(matrixProfile), matrixProfile, 'b');
            hold(gui.profileAx, 'off');
        end
        
        function plotMotifs(gui, motifIdxs)
            if isempty(motifIdxs)
                error('empty motif plot');
            end
            
            % plot top motif pair on data
            gui.plotData();
            
            % checking that the first pair of motifs exists
            if any(isnan(motifIdxs(1 : 2, 1)))
                warning('no valid motifs were provided'); 
                set(gui.dataText, 'string', sprintf('The input time series. No motifs were found.'));
                return;
            end
            
            hold(gui.dataAx, 'on');
            for j = 1:2
                plot(gui.dataAx, motifIdxs(j, 1) : motifIdxs(j, 1) + gui.subLen - 1, gui.data(motifIdxs(j, 1) : motifIdxs(j, 1) + gui.subLen - 1), gui.motifColor(j));
            end
            hold(gui.dataAx, 'off');
            
            % plot motif's neighbors on motif axis
            lims = zeros(3, 2);
            for j = 1:3
                hold(gui.motifAx(j), 'on');
                for k = 3 : size(motifIdxs, 1)
                    if isnan(motifIdxs(k, j))
                        break;
                    end
                    mot = gui.normalizedSubseq(motifIdxs(k, j));
                    plot(gui.motifAx(j), mot, 'color', [0.5 0.5 0.5]); % removed linewidth 2, it was weird having neighbors use a heavier line than motifs
                    lims(j,:) = [min(lims(j, 1), min(mot)); max(lims(j, 2), max(mot))];
                end
                hold(gui.motifAx(j), 'off');
            end
            
            % plot motifs on motif axis
            for j = 1:3
                if any(isnan(motifIdxs(1 : 2, j)))
                    warning('Only %d valid motif set(s) were provided. This is somewhat common with very short time series or long window lengths.', j - 1);
                    break;
                end
                hold(gui.motifAx(j), 'on');
                set(gui.motifText(j), 'string', sprintf(gui.txtTemp{j}, motifIdxs(1, j), motifIdxs(2, j)));
                for k = 1:2
                    mot = gui.normalizedSubseq(motifIdxs(k, j));
                    lims(j,:) = [min(lims(j, 1), min(mot)); max(lims(j, 2), max(mot))];
                    plot(gui.motifAx(j), 1 : gui.subLen, mot, gui.motifColor(k));
                end
                if lims(j, 1) ~= lims(j, 2)
                    ylim(gui.motifAx(j), 1.05 * lims(j,:));
                end
                hold(gui.motifAx(j), 'off');
            end
            set(gui.dataText, 'string', sprintf(['The input time series: ', 'The motifs are color coded (see third panel)']));
        end
        
        function plotDiscords(gui, discordIdxs)
            disclims = zeros(2, 1);
            hold(gui.discordAx, 'on');
            for j = 1:3
                if isnan(discordIdxs(j))
                    warning('Only %d valid discords were provided. This sometimes happens with very short time series or very long window lengths', j - 1);
                    break;
                end
                disc = gui.normalizedSubseq(discordIdxs(j));
                disclims = [min(disclims(1), min(disc)); max(disclims(2), max(disc))];
                plot(gui.discordAx, disc, gui.discordColor(j));
                set(gui.discordText, 'string', sprintf(['The top three discords ', '%d(blue), %d(red), %d(green)'], discordIdxs(1), discordIdxs(2), discordIdxs(3)));
            end
            if disclims(1) ~= disclims(2)
                ylim(gui.discordAx, 1.05 * disclims);
            end
            hold(gui.discordAx, 'off');
        end
        
        function stopCallback(gui)
            gui.shouldHalt_ = true;
            drawnow;
        end
        
        function mainResize(gui)
            figPosition = get(gui.fig, 'position');
            axGap = 38;
            % The zero guard is needed to avoid setting height or position to negative values
            % which will otherwise throw a runtime error.
            axesHeight = max(0, round((figPosition(4) - axGap * 5 - 60) / 6));
            ax_pos = max(figPosition(3) - 60, 0);
            disctxt_pos = max(figPosition(3) - 160, 0);
            set(gui.dataAx, 'position', [30, 5 * axesHeight+5 * axGap + 30, ax_pos, axesHeight]);
            set(gui.profileAx, 'position', [30, 4 * axesHeight+4 * axGap + 30, ax_pos, axesHeight]);
            set(gui.discordAx, 'position', [30, 30, disctxt_pos, axesHeight]);
            set(gui.dataText, 'position',  [30, 6 * axesHeight + 5 * axGap + 30, ax_pos, 18]);
            set(gui.profileText, 'position', [30, 5 * axesHeight + 4 * axGap + 30, ax_pos, 18]);
            set(gui.discordText, 'position', [30, 1 * axesHeight + 30, disctxt_pos, 18]);
            
            % A larger number here places an object higher on the figure.
            for i = 1:3
                set(gui.motifAx(i), 'position', [30, (4 - i) * axesHeight + (4 - i) * axGap + 30, disctxt_pos, axesHeight]);
            end
            for i = 1:3
                set(gui.motifText(i), 'position', [30, (5 - i) * axesHeight + (4 - i) * axGap + 30, disctxt_pos, 18]);
            end
            set(gui.stopBtn, 'position', [ax_pos - 50, 30, 90, 20]);
            drawnow;
        end
        
    end
    
    methods(Access = public)
        function close(gui)
            if isvalid(gui.fig)
                close(gui.fig);
            end
        end
        
        function sh = shouldHalt(gui)
            sh = gui.shouldHalt_;
        end
    end
    
    methods(Static)
        function gui = launchGui(data, mu, invnorm, matrixProfile, motifIdxs, discIdxs, subLen)
            % This is intended as a 1 button catch all to set up everything
            % in a single pass.
            gui = mpgui(data, mu, invnorm, subLen);
            gui.plotProfile(matrixProfile);
            gui.plotMotifs(motifIdxs);
            gui.plotDiscords(discIdxs);
            gui.fig.Visible = 'on';
        end
    end
end