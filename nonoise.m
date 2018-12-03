clear all
close all

load mt_store.csv

load notes.csv

load lattice_store.csv

burn=notes(1);

test=notes(2);

time=notes(3);

len=notes(4);



disp(['burn is ',num2str(burn)])
disp(['test is ',num2str(test)])
disp(['side is ',num2str(len)])
disp(['time is ',num2str(time)])
time_series=zeros(len,len,numel(lattice_store)/len^2);

for b=1:numel(lattice_store)/len^2
for a=1:len
time_series(a,:,b)=lattice_store((b-1)*(len*len)+1+(len*(a-1)):(b-1)*(len*len)+len+len*(a-1));
end
end
over=0;
figure(1)
for b=1:numel(lattice_store)/len^2
   if mod(b,2)==0
   pcolor(time_series(:,:,b))
   end
   xlim([0 len])
   ylim([0 len])
   zlim([0 1])
    pause(0.3)
    if mod(b,50)==0
        prompt = '1 or 2?';
    over = input(prompt);
    end
    if over==1
        break
    end
end

figure(2)
plot(mt_store)
figure(3)
pcolor(time_series(:,:,b))
