clear, close all
processed_folder = 'E:\Dropbox\Research\Traffic\freeway_data\210E\data';
vds = 717635;

% Days in which data set is obtained
dayFirst = datenum('01-Oct-2014');
dayLast = datenum('31-Dec-2014');
days = dayFirst:dayLast;

% Construct PeMS data structure
pemsData = PeMS5minData;
pemsData.load(processed_folder,vds,days);
X = pemsData.get_data_batch_aggregate(vds,days,'var',{'time','flw','dty','spd','occ'},'smooth',false,'fill',false);

% For each day of a week, cluster the days into good and bad days and
% remove bad ones
[X,days]=removeBadDays(X,days);

%% Set parameters
nPointsDay = size(X.flw,1);     % Number of data points
nDays = length(days);           % Number of days
[DayNumber,DayName] = weekday(days);

%%
x = X.flw;
xcurrent = x(1:end-1,:);
xnext = x(2:end,:);
fCurrent = xcurrent;
fNext = xnext;

x = X.occ;
xcurrent = x(1:end-1,:);
xnext = x(2:end,:);
oCurrent = xcurrent;
oNext = xnext;

x = X.spd;
xcurrent = x(1:end-1,:);
xnext = x(2:end,:);
sCurrent = xcurrent;
sNext = xnext;

x = X.time;
x = x - repmat(x(1,1,:),nPointsDay,1,1);
xcurrent = x(1:end-1,:);
xnext = x(2:end,:);
tCurrent = xcurrent;
tNext = xnext;

d = repmat(DayNumber, nPointsDay-1, 1);
d = d(:);

normm = @(x) (x-repmat(mean(x),size(x,1),1))./repmat(std(x),size(x,1),1);

in = [fCurrent, oCurrent, sCurrent, d, tCurrent]';
in = normm(in')';
in(4,:) = in(4,:)*1e1;
out = [fNext, oNext, sNext]';
outM = mean(out,2);
outS = std(out,2);
out = normm(out')';

        
%%
tic
net = feedforwardnet(50,'trainlm');
% net = feedforwardnet([50],'trainscg');
net.trainParam.max_fail = 100;
% net.trainParam.epochs = 00500;
% net.trainParam.Mu = 0.5; 
% net.trainParam.mu_dec = 0.99;
% net.layers{2}.transferFcn = 'tansig'
% net.layers{3}.transferFcn = 'tansig'

net = train(net,in,out,'useparallel','no');
% net = train(net,in,out,'useparallel','no');
toc
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
    for k = 2:nPointsDay
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
col = jet(7);
dtmp = [1,3,7];

figure(1);
hold on;
flw = reshape(X.flw, nt, []);
H=[];
timetmp = linspace(0,24,nt);
for i = dtmp
    idx = find(DayNumber == i);
    idx = idx(5);
    h = plot(timetmp', flw(:,idx),'linewidth',2,'color',col(i,:));
    H = [H;h(1)];
end
[~,tmp]=weekday(dtmp);
legend(H, tmp);
xlabel('Time');
ylabel('Flow (Veh/Hr)');
grid on;
xlim([0,24])
tickValues = 0:1:24;
set(gca,'XTick',tickValues)

figure(2);
hold on;
occ = reshape(X.occ, nt, []);
H=[];
timetmp = linspace(0,24,nt);
for i = dtmp
    idx = find(DayNumber == i);
    idx = idx(5);
    h = plot(timetmp', occ(:,idx),'linewidth',2,'color',col(i,:));
    H = [H;h(1)];
end
[~,tmp]=weekday(dtmp);
legend(H, tmp);
xlabel('Time');
ylabel('Occupancy (%)');
grid on;
xlim([0,24])
tickValues = 0:1:24;
set(gca,'XTick',tickValues)

figure(3);
hold on;
spd = reshape(X.spd, nt, []);
H=[];
timetmp = linspace(0,24,nt);
for i = dtmp
    idx = find(DayNumber == i);
    idx = idx(5);
    h = plot(timetmp', spd(:,idx),'linewidth',2,'color',col(i,:));
    H = [H;h(1)];
end
[~,tmp]=weekday(dtmp);
legend(H, tmp);
xlabel('Time');
ylabel('Speed (mph)');
grid on;
xlim([0,24])
tickValues = 0:1:24;
set(gca,'XTick',tickValues)




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