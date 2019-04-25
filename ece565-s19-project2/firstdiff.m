function [diffcode] = firstdiff(fcc, CONN)
%fcc = [0 2 4 6]
%CONN = 8;
if CONN == 4
    idx = [1 2 3 4 1 2 3 4];
elseif CONN == 8
    idx = [1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8];
end
length = size(fcc);
diffcode = zeros(length(2));
diff = zeros(length(2));
for i = 1:length(2)
    if i == length(2)
        cur_val = fcc(1,i);
        nex_val = fcc(1,1);
    else
        cur_val = fcc(1,i);
        nex_val = fcc(1,i+1);
    end
    count = 0;
    diff = 0;
    for j = 1:CONN
        if idx(cur_val+1 + j - 1) ~= idx(nex_val+1)
            count = count + 1;
        elseif idx(cur_val+1 + j - 1) == idx(nex_val+1)
            diff(i) = count;
        end
    end  
end
diffcode = diff;
end