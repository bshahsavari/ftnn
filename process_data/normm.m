function x = normm(x) 
% Normalize a matrix to have zero mean columns with unity standard
% deviation
x=(x-repmat(mean(x),size(x,1),1))./repmat(std(x),size(x,1),1);