function [group1 group2] = make_group12(matchList)


nMatch = size(matchList,1);

featList = matchList(:,1);
featList_unique = unique(featList);
nGroup = length(featList_unique);
group = logical(sparse(nMatch,nGroup));
for i=1:nGroup
    group( find(featList == featList_unique(i)), i) = true;
end
group1 = group;

featList = matchList(:,2);
featList_unique = unique(featList);
nGroup = length(featList_unique);
group = logical(sparse(nMatch,nGroup));
for i=1:nGroup
    group( find(featList == featList_unique(i)), i) = true;
end
group2 = group;