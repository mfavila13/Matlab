function [diffcode] = firstdiff(fcc, CONN)

if CONN == 4
    idx = [0 1 2 3 0 1 2 3];
elseif CONN == 8
    idx = [0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7];
end
length = size(fcc);
diffcode = zeros([1 length(2)]);
diff = zeros([1 length(2)]);
% set indexing
for i = 1:length(2)
    if i == length(2)
        cur_val = fcc(i);
        nex_val = fcc(1);
    else
        cur_val = fcc(i);
        nex_val = fcc(i+1);
    end
    count = 0;
    %count number difference
    for j = 1:CONN
        if idx(cur_val+1 + j - 1) ~= idx(nex_val+1)
            count = count + 1;
        elseif idx(cur_val+1 + j - 1) == idx(nex_val+1)
            diff(i) = count;
        end
    end  
end
diffcode = diff
end