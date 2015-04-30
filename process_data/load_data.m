processed_folder = 'E:\Dropbox\Research\Traffic\freeway_data\210E\data';
X.good_days = get_good_days();
vds = [717635];
days = X.good_days;
dayFirst = datenum('01-Oct-2014');
dayLast = datenum('31-Dec-2014');
days = dayFirst:dayLast;

pemsData = PeMS5minData;
pemsData.load(processed_folder,vds,days)

X=pemsData.get_data_batch_aggregate(vds,days,'var',{'time','flw','dty','spd','occ'},'smooth',false,'fill',false)
%%
nd = length(days);
figure; 
xmin = min(min(X.dty));
xmax = max(max(X.dty));
ymin = min(min(X.flw));
ymax = max(max(X.flw));
for i = 1:nd
    plot(X.dty(:,1,i),X.flw(:,1,i),'.');
    grid on
    axis([xmin,xmax,ymin,ymax]);
    pause
end
    
%%
nt = size(X.flw,1);
x = X.flw;
xcurrent = x(1:end-1,:);
xcurrent = xcurrent(:);
xnext = x(2:end,:);
xnext = xnext(:);
fCurrent = xcurrent;
fNext = xnext;


x = X.occ;
xcurrent = x(1:end-1,:);
xcurrent = xcurrent(:);
xnext = x(2:end,:);
xnext = xnext(:);
oCurrent = xcurrent;
oNext = xnext;

x = X.spd;
xcurrent = x(1:end-1,:);
xcurrent = xcurrent(:);
xnext = x(2:end,:);
xnext = xnext(:);
sCurrent = xcurrent;
sNext = xnext;

x = X.time;
x = x - repmat(x(1,1,:),nt,1,1);
xcurrent = x(1:end-1,:);
xcurrent = xcurrent(:);
xnext = x(2:end,:);
xnext = xnext(:);
tCurrent = xcurrent;
tNext = xnext;

[DayNumber,DayName] = weekday(days);
d = repmat(DayNumber, nt-1, 1);
d = d(:);

normm = @(x) (x-repmat(mean(x),size(x,1),1))./repmat(std(x),size(x,1),1);

in = [fCurrent, oCurrent, sCurrent, d, tCurrent];
% in = [oCurrent, sCurrent, d, tCurrent];
in = normm(in);
out = [fNext, oNext, sNext];
% out = [oNext, sNext];
outM = mean(out);
outS = std(out);
out = normm(out);

%%
clear y
idx = 88*287+1;
y(:,1) = in(idx,1:3)';
for k = 2:nt
    tmp = [y(:,k-1);in(k-1+idx-1,4);in(k-1+idx-1,5)];
    if ~mod(k,50)
        tmp = in(k-1+idx-1,:)';
    end
    y(:,k) = myNeuralNetworkFunction(tmp);
end

y = y.*repmat(outS',1,size(y,2)) + repmat(outM',1,size(y,2));

figure; plot([y(3,:)' sCurrent(idx:idx+287,1)],'.-');
figure; plot([y(2,:)' oCurrent(idx:idx+287,1)],'.-');
figure; plot([y(1,:)' fCurrent(idx:idx+287,1)],'.-');
% y(:,1) = in(1,1:2)';
% for k = 2:nt
%     tmp = [y(:,k-1);in(k-1,3);in(k-1,4)];
%     if ~mod(k,1)
%         tmp = in(k-1,:)';
%     end
%     y(:,k) = myNeuralNetworkFunction4(tmp);
% end
% 
% y = y.*repmat(outS',1,size(y,2)) + repmat(outM',1,size(y,2));
% 
% figure; plot([y(2,:)' sCurrent(1:288,1)],'.-');
% figure; plot([y(1,:)' oCurrent(1:288,1)],'.-');
%%


%%
clear y
y(:,1) = in2(1,1:3)';
yold = y;
m=3
j=1
for k = 2:m:nt
    tmp = [yold;d(k-1);m];
    if ~mod(k,m)
        tmp = [in(k,:)';m];
    end
    y(:,j) = myNeuralNetworkFunction3(tmp);
    j = j+1;
    yold = y(:,end);
end

y = y.*repmat(outS',1,size(y,2)) + repmat(outM',1,size(y,2));

figure; plot([y(3,:)' sNext(1:m:288,1)],'.-');