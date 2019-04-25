function [minimum] = min_mag(fcc)

length = size(fcc);
all_fcc = fcc;
for j = 2:length(2)
    all_fcc(j,1:length(2)-1) = all_fcc(j-1,2:length(2));
    all_fcc(j,length(2)) = all_fcc(j-1,1);
end
index = 1;
% find minimum of each column, if not minimum, set to all 9's
for i = 1:size(all_fcc,2)
    min_val = min(all_fcc(:,i));
    for k = 1:size(all_fcc,2)
        if all_fcc(k,i) ~= min_val
            all_fcc(k,:) = 9*ones(1,size(all_fcc,2));
        end
    end
end
% minimum is row that is not all 9's
for q = 1:size(all_fcc,2)
    if all_fcc(q,1) ~= 9
        minimum = all_fcc(q,:)
    end
end
%minimum = all_fcc
end