function [minimum] = min_mag(fcc)
%{
min_int = str2double(sprintf('%d',round(fcc)));

shift = zeros(length)
for i = 1:length-1
    shift(i) = fcc(i+1)
end
shift(length) = fcc(1)
integer = str2double(sprintf('%d',round(a)));
if integer < min_int
    minimum = shift
%}
fcc
length = size(fcc);
all_fcc = fcc;
for j = 2:length(2)
    all_fcc(j,1:length(2)-1) = all_fcc(j-1,2:length(2));%fcc(2:length(2))
    all_fcc(j,length(2)) = all_fcc(j-1,1);%fcc(1);
end
index = 1;
for i = 1:size(all_fcc,2)
    min_val = min(all_fcc(:,i));
    for k = 1:size(all_fcc,2)
        if all_fcc(k,i) ~= min_val
            all_fcc(k,:) = 9*ones(1,size(all_fcc,2));
        end
    end
end
for q = 1:size(all_fcc,2)
    if all_fcc(q,1) ~= 9
        minimum = all_fcc(q,:);
    end
end
%minimum = all_fcc
end