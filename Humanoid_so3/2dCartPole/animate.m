function animate(z, sim_dt,p)
    f3 = figure;
    figure(f3);
    axis equal
    frames = length(z);
    cart_w = 0.4;
    cart_h = 0.2;

    for i = 1:frames
        if mod(i, 10)==0
            positions = keypoints_position(z(:,i), p);
            x = positions(1,:);
            y = positions(2,:);
            figure(f3);
            plot(x,y, 'LineWidth',5);
            rectangle('position', [x(2)-0.4+cart_w/2 -cart_h/2 cart_w cart_h])
            %plotcube([0.4 0.4 0.2],[x(2)-0.2 -0.2 -0.2],.5,[1 0 0]);
            hold on
            circle(x(1),y(1), 0.05)
            axis([-5 5 -5 5]);
            pause(0.0002);
            hold off
        
        end


    end

end

function circle(x,y,r)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp, 'LineWidth',5);
end
