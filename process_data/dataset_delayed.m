%%out = dataset_delayed(data,N,dmin,dmax)
% Crates normalized input and output data for NN.
% data is a matrix with columns: flow, occ, spd, day, time
% N number of rows of input and output
% dmin: minimum number of delay steps between input and output of NN
% dmax: maximum   "    "   "    ...
function out = dataset_delayed(data,N,dmin,dmax)

d = data;
% N = 26e3;1e5;
% dmin = 1;6;
% dmax = 1;17;

nsamples = size(d,1);
input = [];
output = [];
numData4EachDelay = floor(N/(dmax-dmin +1));

% Indices of data: flow, occ, speed, day, time
fi      = 1;
oi      = 2;
si      = 3;
dayi    = 4;
ti      = 5;

% Construct input/output data with N rows
while size(output,1) < N
    for delay = dmin:dmax
        
        % Choose random data instances
        idx = randperm(nsamples-delay);
        % For the last delay, draw as many as needed to make the total number of rows equal to N 
        if delay == dmax
            numData4EachDelay = min(length(idx),N - size(input,1));
        end
        idx = idx(1:numData4EachDelay);

        % Remove negative delays
        delayVector = d(idx+delay,ti)-d(idx,ti);
        negIdx = find(delayVector < 0);
        idx(negIdx) = [];
        delayVector(negIdx) = [];

        % Construct input training data: flow, occ, speed, day, time, delay
        tmp = [d(idx,[fi,oi,si,dayi,ti]) delayVector];
        input = [input; tmp];

        % Index of delayed data
        idx = idx+delay;

        % Construct output training data
        tmp = [d(idx,[fi,oi,si])];
        output = [output;tmp];
    end
end

% Shuffle data
shufIdx = randperm(N);
input = input(shufIdx,:);
output = output(shufIdx,:);

% Remove any extra row
input = input(1:N,:);
output = output(1:N,:);

% Normalize data
inMean = mean(input);
inStd = std(input);
input = normm(input);
outMean = mean(output);
outStd = std(output);
output = normm(output);

% Construct output argument
out = struct('output',output, 'input', input, 'inMean', inMean, 'inStd',inStd,'outMean', outMean, ...
    'outStd', outStd);