net = feedforwardnet(10,'trainlm');
net.trainParam.max_fail = 100;
net = train(net,in',out');
view(net)
y = net(x);
perf = perform(net,y,t)