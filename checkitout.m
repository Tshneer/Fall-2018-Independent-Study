clear all;

close all

% 
load noise_range.csv
%load mean_pop.csv
%load mean_square.csv
load lattice_store.csv
load sus1.csv
load sus2.csv

load order.csv

load lattice.csv

load notes.csv

len=notes(5);
lat2d=zeros(len,len,numel(noise_range)); 

%sus=sus1-sus2.*sus2;

sus=sus1-sus2.*sus2;



for b=1:numel(noise_range)
for a=1:len
lat2d(a,:,b)=lattice_store((b-1)*(len*len)+1+(len*(a-1)):(b-1)*(len*len)+len+len*(a-1));
end
end
ord = ['The order is ', num2str(notes(1)),'.'];

global_range=notes(4);

%loc = ['The local noise strength is ',num2str(notes(2)),'.'];
%glob = ['The global noise strength is ',num2str(notes(3)),'.'];
burn = ['The burn is ',num2str(notes(1)),'.'];
test = ['The measure time is ',num2str(notes(2)),'.'];
run = ['The run time is ',num2str(notes(3)),' seconds.'];


disp(burn)
disp(test)
disp(run)
%figure(1)
%pcolor(lat2d)

% figure(2)
% plot(order,'*')
% 
% std=sqrt( mean_square - mean_pop.^2);





reduced_noise=2.85.*len.*(noise_range-0.05074)./noise_range;



% plot(reduced_noise,order,'.')
% 
% xlabel('Local Noise Level')
% ylabel('Order Level')
% xlim([-4,4])
% subplot(3,1,2)
% plot(order,mean_pop./((notes(2))*(len^2)),'.')
% xlabel('Order Level')
% ylabel('Time average of Population Density')

sus=len^2*sus.*52.243984.*len.^(-7/4);
row_length=numel(sus)/global_range;
sus_mat=zeros(global_range,row_length);
for row=1:(global_range)
sus_mat(row,:)=sus((row-1)*row_length+1:row*row_length);
end
order_mat=zeros(global_range,row_length);
for row=1:(global_range)
order_mat(row,:)=order((row-1)*row_length+1:row*row_length);
end


for yay=1:global_range
    figure(yay)
 subplot(2,1,1)

plot(reduced_noise,order_mat(yay,:),'.')
 subplot(2,1,2)
 plot(reduced_noise,sus_mat(yay,:),'.')
 title(['Global Strength is: ',num2str((yay-1)*0.01),])
 xlabel('Local Noise Level')
 ylabel('Sus')
 xlim([-4,4])

end

max_sus=zeros(global_range,1);
[M,I] =max(sus_mat');
% for a=1:(global_range+1)
% max_sus(a)=max(sus_mat,[],11);
% end
figure(yay+1)
title('local noise level of max variance ')
ylabel('local noise level')
plot(noise_range(I),'--')


figure(yay+2)

plot(M,'--')
title('max variance')


for b=1:numel(noise_range)
figure(b)
pcolor(lat2d(:,:,b))
end