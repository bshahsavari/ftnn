clear all; close all
processed_folder = 'E:\Dropbox\Research\Traffic\freeway_data\210E\data';
vds = [717635];

% Use this for selecting good days
% days = get_good_days();

% Use this for selecting all days
dayFirst = datenum('01-Oct-2014');
dayLast = datenum('31-Dec-2014');
days = dayFirst:dayLast;

pemsData = PeMS5minData;
pemsData.load(processed_folder,vds,days)
X=pemsData.get_data_batch_aggregate(vds,days,'var',{'time','flw','dty','spd','occ'},'smooth',false,'fill',false)


%% Parse data
nt = size(X.flw,1);

flw = X.flw(:);
occ = X.occ(:);
spd = X.spd(:);

x = X.time;
x = x - repmat(x(1,1,:),nt,1,1);
t = x(:);

[DayNumber,DayName] = weekday(days);
d = repmat(DayNumber, nt, 1);
d = d(:);

% This is the main data container
data = [flw, occ, spd, d, t];

%% 
out = dataset_delayed(data,2.6e4,1,1);
in1 = out.input;
out1 = out.output;
%% Simulate NN
day1 = 2;
m = 1;
f = @n4
K = 280
simulate_nn(data,nt,out,day1,f,m,K)

%%
break
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

%% Plot flow 
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