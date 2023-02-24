function cs = buildcs(num_steps,swingt)
    cs_len = num_steps*5+1;

    cs = ones(4, cs_len);
    flag = 1;
    for i= 2:swingt:cs_len
        if flag ==1
            cs(1,i:i+swingt-1)= [1 1 1 1 0];    %leg1 right heel 
            cs(2,i:i+swingt-1)= [1 1 1 1 1];    %leg1 right toe 
            cs(3,i:i+swingt-1)= [0 0 0 0 1];    %leg2 left heel 
            cs(4,i:i+swingt-1)= [0 0 0 0 0];    %leg2 left toe 
            flag = 0;
        else
            cs(1,i:i+swingt-1)= [0 0 0 0 1];    %leg1 right heel 
            cs(2,i:i+swingt-1)= [0 0 0 0 0];    %leg1 right toe 
            cs(3,i:i+swingt-1)= [1 1 1 1 0];    %leg2 left heel 
            cs(4,i:i+swingt-1)= [1 1 1 1 1];    %leg2 left toe 
            flag = 1;
        end

    end

end
