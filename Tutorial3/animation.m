% enable:       Returns without action if false [true/false]
% t:            Time series [t1; t2; ...; tn]
% tEnd:         Use to cut animation short. Will display until tn if tEnd > tn
% targetFps:    Targed display fps. Rounds to nearest possible value based
%               on delta time and speed
% spd:          Animation speed (x)
% vewBox:       Viewport size [xMin, xMax, yMin, yMax]
% drawCallback: Callback function for drawing information.
%               Signature: drawCallback(figureHandle, axisHandle, drawData, t, headline, fps)
% drawData:     Data passed to drawCallback [data(t1), data(t2), ..., data(tn)]
function animation(enable, t, tEnd, targetFps, spd, headline, viewBox, drawCallback, drawData)
    if ~enable
        return;
    end
    
    dt = t(2) - t(1);
    nSamples = length(t);
    stp = round(spd/(dt*targetFps));
    actDt = stp*dt/spd;
    
    figureHandle = figure;
    xlabel('X [m]')
    ylabel('Y [m]')
    figureHandle.WindowState = 'maximized';
    axisHandle = gca;
    axis equal;
    axis(viewBox);
    hold on;
    % Small delay to make sure the window has been properly created before
    % starting to draw
    pause(2);
    for i = 1:stp:nSamples
        if i > 1
            fps = 1/toc;
        else
            fps = 0;
        end
        tic;

        % Clear axis for next draw
        cla(axisHandle);
        
        % Draw objects
        drawCallback(figureHandle, axisHandle, drawData(:, i), t(i), headline, fps);

        if t(i) >= tEnd
            break;
        end
        drawTime = toc;
        delayTime = dt*stp/spd - drawTime;
        if delayTime < 0
            error(['Could not maintain desired fps (', num2str(1/actDt), '). Drawtime: ', num2str(drawTime), 's'])
        end
        pause(delayTime);
    end
end