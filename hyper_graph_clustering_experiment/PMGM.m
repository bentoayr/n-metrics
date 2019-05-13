function [X] = PMGM(rawMat,nodeCnt,graphCnt)
n = nodeCnt;

N = max([nodeCnt,graphCnt]);
cutoffeig_num = n;

[u,~,~] = svd(rawMat);
cu = u(:,1:cutoffeig_num);

B={};
for i = 1:n:N*n
    B{(i+(n-1))/n} = cu(i:(i+n)-1,:);
end

recPM=cell(N,1);
recDiscMat = zeros(nodeCnt*graphCnt,nodeCnt*graphCnt);
I = eye(n);
for i = 1:N
    tmp =  B{i}*B{1}';
    [tmpU,~,tmpV]=svd(tmp);
    recPM{i} = tmpU*tmpV';

    [c,~] = munkres(-recPM{i}');
    recPM{i} = I(:,c);

end
for i = 1:N
    for j = i+1:N
        recDiscMat((i-1)*n+1:i*n,(j-1)*n+1:j*n) = recPM{i}*recPM{j}';       
    end
end
X = recDiscMat + recDiscMat' + eye(nodeCnt*graphCnt);