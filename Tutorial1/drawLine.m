function drawLine(axHandle, pos, rot, length, color, thickness)
%    rot = [cos(angle); sin(angle)];
    offs = length/2*rot;
    
    pt0 = pos - offs;
    pt1 = pos + offs;
    
    plot(axHandle, [pt0(1), pt1(1)], [pt0(2), pt1(2)] , 'LineWidth', thickness, 'Color', color);
end