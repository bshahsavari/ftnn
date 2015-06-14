function [X,days]=removeBadDays(X,days)

DayNumber = weekday(days);
nPointsDay = size(X.flw,1);     % Number of data points

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
    hr(:,1)=plot(reshape(x(:,1,goodDays),nPointsDay,[]),'b');
    hold on
    hr = [hr; plot(reshape(x(:,1,badDays),nPointsDay,[]),'r')];
    legend(hr([1 end]),'Good days', 'Bad days')

    X.flw(:,:,badDays)=[];
    X.dty(:,:,badDays)=[];
    X.spd(:,:,badDays)=[];
    X.occ(:,:,badDays)=[];
    X.time(:,:,badDays)=[];
    days(badDays)=[];
    DayNumber = weekday(days);
%     pause
end