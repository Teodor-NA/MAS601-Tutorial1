function drawBox(axHandle, pos, angle, height, width, color, thickness)
    ca = cos(angle);
    sa = sin(angle);
    rot = [ca, -sa; sa, ca];
    xOffs = rot*[width/2; 0];
    yOffs = rot*[0; height/2];
    
%     rotX = [cos(angle); sin(angle)];
%     rotY = [cos(angle + pi/2); sin(angle + pi/2)];
%     xOffs = width/2*rotX;
%     yOffs = height/2*rotY;
     
    drawLine(axHandle, pos - xOffs, angle + pi/2, height, color, thickness); 
    drawLine(axHandle, pos + xOffs, angle + pi/2, height, color, thickness);
    drawLine(axHandle, pos - yOffs, angle, width, color, thickness); 
    drawLine(axHandle, pos + yOffs, angle, width, color, thickness);
end