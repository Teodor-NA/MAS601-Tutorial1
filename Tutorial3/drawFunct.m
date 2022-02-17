function drawFunct(figureHandle, axisHandle, drawData, t, headline, fps)
global L1 L2 L3 L4
    title(axisHandle, ...
        [headline, ' --- t: ', num2str(t, '%.1f'), '[s] --- ', num2str(fps, '%.1f'), ' fps']);

    drawBox(axisHandle, [0; 0], 0, 0.2, 0.2, 'k', 4);
    drawBox(axisHandle, [L1; 0], 0, 0.2, 0.2, 'k', 4);
    drawBox(axisHandle, drawData(1:2), drawData(3), 0.2, L2, 'r', 4);
    drawBox(axisHandle, drawData(4:5), drawData(6), 0.2, L3, 'b', 4);
    drawBox(axisHandle, drawData(7:8), drawData(9), 0.2, L4, 'g', 4);
end