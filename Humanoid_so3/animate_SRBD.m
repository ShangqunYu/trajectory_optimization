function animate_SRBD(xfb, rot, cfootheel, cfoottoe, cs, r, P_des, R_des)
    f3 = figure;
    hold off
    figure(f3);
    axis equal
    %axis([-1 1 -1 1 -1 1]);
    frames = length(xfb);

    for i = 1:frames
        
        figure(f3);
        clf
        grid on
        % plotting the floor
        patch([3,-3,-3,3],[3,3,-3,-3],[0,0,0,0],[0.5 0.5 0.5], 'FaceAlpha',.3);

        plotSRBD([r r r],xfb(1:3,i)', rot(:,:,i));
        plotSRBD([r r r],P_des(:,i)', R_des(:,:,i), 0.2, [0 0.5 0]  );
        hold on
        if i < frames
            if cs(1,i)==1
                plot3(cfootheel(1,i),cfootheel(2,i),cfootheel(3,i), 'o','Color','r');  % right heel "o"
            end
            if cs(2,i)==1
                plot3(cfoottoe(1,i),cfoottoe(2,i),cfoottoe(3,i), "^",'Color','r');  % right toe "^"
            end
            if cs(3,i)==1
                plot3(cfootheel(4,i),cfootheel(5,i),cfootheel(6,i), 'o','Color','b');  % left heel
            end
            if cs(4,i)==1
                plot3(cfoottoe(4,i),cfoottoe(5,i),cfoottoe(6,i), "^",'Color','b');  % left toe "^"
            end
        end
        axis([-1 4 -3 3 0 2]);
        pause(0.5);
        disp(i);
    end



    



end
