function [ idx ID ] = make_groups( group )

nGroup = size(group,2);
idx = [];
ID = [];
for i = 1:nGroup
    tmp = find(group(:,i));
    idx = [idx; tmp];
    ID = [ID; i*ones(size(tmp))];
end