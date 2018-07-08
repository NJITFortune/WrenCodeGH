function idx = wreK(struct, knum)

cols(1,:)='r*'; cols(2,:)='b*'; cols(3,:)='m*'; cols(4,:)='g*'; 
cols(5,:)='c*'; cols(6,:)='k*'; cols(7,:)='y*'; cols(8,:)='ro';

%% First things first - pick off the syllables with the highest frequencies




a(:,1) = [struct.sylen];
a(:,2) = [struct.trace_peakf];
a(:,3) = [struct.trace_peakt];

%figure(2); clf; plot3(a(:,1), a(:,2), a(:,3), 'k*');
figure(2); clf; plot(a(:,1), a(:,2), 'k*');

idx = kmeans(a, knum);

hold on;
for i=1:knum;
    b = idx == i;
%    figure(2); plot3(a(b,1), a(b,2), a(b,3), cols(i,:));
    figure(2); plot(a(b,1), a(b,2), cols(i,:));
end;



