function simulate_nn(data,nt,train,daystart,f,m,K)

fi      = 1;
oi      = 2;
si      = 3;
dayi    = 4;
ti      = 5;
deli    = 6;

dayend = daystart+1;
idx = (daystart-1)*(nt-1)+1;
idxend = (dayend-1)*(nt-1);
input = data(idx:m:idxend,:);
inputs = (input - repmat(train.inMean(1:5),size(input,1),1))...
    ./repmat(train.inStd(1:5),size(input,1),1);

% delM = train.inMean(deli);
% delS = train.inStd(deli);

y(:,1) = inputs(1,[fi oi si])';
min5 = 1/24/60*5;
for k = 2:size(input,1)
    tmp = [y(:,k-1);inputs(k-1,[dayi ti])'];
    if ~mod(k,K)
        tmp = [inputs(k-1,:)'];
    end
    y(:,k) = f(tmp);
end

y = y.*repmat(train.outStd',1,size(y,2)) + repmat(train.outMean',1,size(y,2));

tim = linspace(0,24,size(y,2));
figure; plot(tim,[y(si,:)' input(:,si)],'.-');
figure; plot(tim,[y(fi,:)' input(:,fi)],'.-');
figure; plot(tim,[y(oi,:)' input(:,oi)],'.-');