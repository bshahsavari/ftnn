clear all; close all;
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


[DayNumber,DayName] = weekday(days);
nt = size(X.flw,1);

%%
nd = length(days);
% figure; 
% xmin = min(min(X.dty));
% xmax = max(max(X.dty));
% ymin = min(min(X.flw));
% ymax = max(max(X.flw));
% for i = 1:nd
%     plot(X.dty(:,1,i),X.flw(:,1,i),'.');
%     grid on
%     axis([xmin,xmax,ymin,ymax]);
%     pause
% end
    
%%
for daytmp = 1:7
    fos = [ (X.flw - mean(X.flw(:)))/std(X.flw(:)); 
            (X.occ - mean(X.occ(:)))/std(X.occ(:));
            (X.spd - mean(X.spd(:)))/std(X.spd(:))];
    fos = reshape(fos,size(fos,1),[]);
    didx = find(DayNumber == daytmp);
    [tmp,~,sumd] = kmeans(fos(:,didx)',2,'Replicates',10);
    goodDays = didx(tmp==mode(tmp));
    badDays = didx(tmp~=mode(tmp));
    sum(sumd)

    x = X.spd;
    clear hr
    figure;
    hr(:,1)=plot(reshape(x(:,1,goodDays),nt,[]),'b');
    hold on
    hr = [hr; plot(reshape(x(:,1,badDays),nt,[]),'r')];
    legend(hr([1 end]),'Good days', 'Bad days')

    X.flw(:,:,badDays)=[];
    X.dty(:,:,badDays)=[];
    X.spd(:,:,badDays)=[];
    X.occ(:,:,badDays)=[];
    X.time(:,:,badDays)=[];
    days(badDays)=[];
    [DayNumber,DayName] = weekday(days);

%     pause
end


%%
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

d = repmat(DayNumber, nt-1, 1);
d = d(:);

normm = @(x) (x-repmat(mean(x),size(x,1),1))./repmat(std(x),size(x,1),1);

in = [fCurrent, oCurrent, sCurrent, d, tCurrent]';
% in = [oCurrent, sCurrent, d, tCurrent];
in = normm(in')';
in(4,:) = in(4,:)*1e1;
out = [fNext, oNext, sNext]';
% out = [oNext, sNext];
outM = mean(out');
outS = std(out');
out = normm(out')';





        
%%
% didx = find(d==7);
% x = oCurrent;
% x = x(didx,1);
% x = reshape(x, nt-1, []);
% idx = kmeans(x', 2);
% idx1 = find(idx == mode(idx));
% 
% x = sCurrent;
% x = x(didx,1);
% x = reshape(x, nt-1, []);
% idx = kmeans(x', 2);
% idx2 = find(idx == mode(idx));
% 
% x = fCurrent;
% x = x(didx,1);
% x = reshape(x, nt-1, []);
% idx = kmeans(x', 2);
% idx3 = find(idx == mode(idx))



%%
net = feedforwardnet([40,6],'trainlm');
% net = feedforwardnet([50],'trainscg');
net.trainParam.max_fail = 2000;
net.trainParam.epochs = 250;
% net.layers{2}.transferFcn = 'tansig'
% net.layers{3}.transferFcn = 'tansig'

net = train(net,in,out,'useparallel','yes');
% net = train(net,in,out,'useparallel','no');

view(net)

y = net(in);
perf = perform(net,y,out)

%%

dn =2;
DayNumber(dn)
dtmp = find(DayNumber ==DayNumber(dn),5)
dtmp = [1:3]
for dn=dtmp
    clear y
    dn
    idx = (dn-1)*287+1;
    y(:,1) = in(1:3,idx);
    for k = 2:nt
        tmp = [y(:,k-1);in(4,k-1+idx-1);in(5,k-1+idx-1)];
        if ~mod(k,288)
            tmp = in(:,k-1+idx-1);
        end
        y(:,k) = net(tmp);
    end

    y = y.*repmat(outS',1,size(y,2)) + repmat(outM',1,size(y,2));

    figure(1); hold on; plot([y(3,:)' sCurrent(idx:idx+287,1)],'.-');
    legend('yh','y')
    figure(2); hold on; plot([y(2,:)' oCurrent(idx:idx+287,1)],'.-');
    figure(3); hold on; plot([y(1,:)' fCurrent(idx:idx+287,1)],'.-');
    drawnow
end

break;
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
dataset_delayed
inM2 = mean(in2);
inS2 = std(in2);
in2 = normm(in2);
outM2 = mean(out2);
outS2 = std(out2);
out2 = normm(out2);

delS = inS2(6);
delM = inM2(6);

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