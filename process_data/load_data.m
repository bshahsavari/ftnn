processed_folder = 'E:\Dropbox\Research\Traffic\freeway_data\210E\data';
X.good_days = get_good_days();
vds = [717635];
days = X.good_days;

pemsData = PeMS5minData;
pemsData.load(processed_folder,vds,days)

X=pemsData.get_data_batch_aggregate(vds,days,'var',{'time','flw','dty','spd','occ'},'smooth',false,'fill',false)
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
% 
% x = X.time;
% xcurrent = x(1:end-1,:);
% xcurrent = xcurrent(:);
% xnext = x(2:end,:);
% xnext = xnext(:);
% tCurrent = xcurrent;
% tNext = xnext;

tCurrent = [1:287]';
tCurrent = repmat(tCurrent,36,1);
tNext = tCurrent+1;

[DayNumber,DayName] = weekday(days);
d = repmat(DayNumber, 287, 1);
d = d(:);

normm = @(x) (x-repmat(mean(x),size(x,1),1))./repmat(std(x),size(x,1),1);

in = [fCurrent, oCurrent, sCurrent, d];
in = normm(in);
out = [fNext, oCurrent, sCurrent];
outM = mean(out);
outS = std(out);
out = normm(out);




