function I = MutualInformation(X,Y)
 
JbinSize = [round(range(X)/std(X)*5),round(range(Y)/std(Y)*5)];
XbinSize = round(range(X)/std(X)*5);
YbinSize = round(range(Y)/std(Y)*5);

    I = Entropy(X, XbinSize) + Entropy(Y, YbinSize) - JointEntropy(X, Y, JbinSize);
 
end
