d = in;
N = 1e5;
dmin = 6;
dmax = 17;

dd = dmax-dmin +1;
nsamples = size(d,1);
in2 = [];
out2 = [];
ndata = floor(N/dd);

% Indices of data
fi      = 1;
oi      = 2;
si      = 3;
dayi    = 4;
ti      = 5;

for delay = dmin:dmax
    % Choose "ndata" random data instances
    idx = randperm(nsamples-delay);
    if delay == dmax
        ndata = N - size(in2,1);
    end
    idx = idx(1:ndata);
    
    % Construct input training data
    delayVector = repmat(delay,length(idx),1);
    tmp = [d(idx,[fi,oi,si,dayi,ti]) delayVector];
    in2 = [in2; tmp];
    
    % Index of delayed data
    idx = idx+delay;
    
    % Construct output training data
    tmp = [d(idx,[fi,oi,si])];
    out2 = [out2;tmp];
end
