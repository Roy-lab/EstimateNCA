function getCorr()

tfas={};
for i=0:9
	a=importdata(sprintf('../tfa/tfa0.010_%d.txt',i));
	tfas=[tfas,a.data];
	tflen = size(a.data,1);
end
allq=[];
for t=1:tflen
	allc=[];
	for i=1:10
		t1=tfas{i};
		for j=i+1:10
			t2=tfas{j};
			c=corr(t1(t,:)',t2(t,:)');
			c=abs(c);
			allc=[allc;c];
		end
	end
	q=quantile(allc,3);
	allq=[allq;q];
end
[~,ids] = sort(allq(:,2));
allq=allq(ids,:);
f=figure('paperunit','inches','paperposition',[0 0 6 4],'visible','off');
plot(allq);
ylabel('Abs Corr');
legend({'25%','50%','75%'},'location','nw');
print(f,'../corr_stability.png','-dpng');
close(f);
