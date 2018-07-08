function H = Entropy(y, binSize)
    % Calculate the entropy for an integer value of y
 
    % I suggest: binSize=range(y)/std(y)*5;
    % Generate the histogram
    [n, xout] = hist(y, binSize);
 
    % Normalize the area of the histogram to make it a pdf
    n = n / sum(n);
    b = xout(2) - xout(1);
 
    % Calculate the entropy
    indices = n ~= 0;
    H = -sum(n(indices) .* log2(n(indices)) .* b);
    
end