function simulate_nn(data,nt,train)

fi      = 1;
oi      = 2;
si      = 3;
dayi    = 4;
ti      = 5;
deli    = 6;

daystart = 1;
dayend = 2;
m = 1;
f = @n3
clear y
idx = (daystart-1)*(nt-1)+1;
idxend = (dayend-1)*(nt-1);
input = data(idx:m:idxend,:);
inputs = (input - repmat(train.inMean(1:5),size(input,1),1))...
    ./repmat(train.inStd(1:5),size(input,1),1);

delM = train.inMean(deli);
delS = train.inStd(deli);

y(:,1) = inputs(1,[fi oi si])';
min5 = 1/24/60*5;
for k = 2:size(input,1)
    tmp = [y(:,k-1);inputs(k-1,[dayi ti])'; (m*min5-delM)/delS];
    if ~mod(k,60)
        tmp = [inputs(k-1,:)'; (m*min5-delM)/delS];
    end
    y(:,k) = f(tmp);
end

y = y.*repmat(train.outStd',1,size(y,2)) + repmat(train.outMean',1,size(y,2));

tim = linspace(0,24,size(y,2));
figure; plot(tim,[y(si,:)' input(:,si)],'.-');
figure; plot(tim,[y(fi,:)' input(:,fi)],'.-');
figure; plot(tim,[y(oi,:)' input(:,oi)],'.-');