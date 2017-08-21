function test_5hz

n=1; test=zeros(11999,34);
for i = 1:34 % each trial
    if exist(int2str(i),'file')==2
        temp=importdata(int2str(i));
        if n~1
            n=n+1;
        end
    end
    test(:,n) = temp(1:11999,1); % I put row 1:points in case  a file was recorded twice in Nguyet
end

duration_ms=11999/20000*1000;
step=duration_ms/11999;
time_axis=0:step:duration_ms-step; time_axis=time_axis';

m=1;
for j = 1:3
    figure;
   for i = 1:9 
      subplot(3,3,i);
      plot(time_axis,test(:,m))
      m=m+1;
   end
end

end