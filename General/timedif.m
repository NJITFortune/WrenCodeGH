function out = timedif

	a = ginput(2);
	a = a(:,1);

	out = abs(a(1) - a(2));
	
