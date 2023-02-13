clc
clear
clf
load("linearity_exercise.mat")
for  i = 1: size(ff_images,3)
     temp = double(ff_images(:,:,i));
     ind = find(temp>0);
     med_val(i) = median(temp(ind));
     avg_val(i) = mean(temp(ind));
     std_dev(i) = std(temp(ind));

end
ph(1) = plot(itime(1:30),med_val(1:30),'b-');
hold on
ph(2) = errorbar(itime(1:30),avg_val(1:30),std_dev(1:30),'Color','red');
xlabel('Integration(ms)');
ylabel('response(DN)');
set(gcf,'Color','w');
grid on
legend([ph(1),ph(2)],'Median response','Average Response')
backstr = [''];
for i = 1:size(ff_images,1)
    for j = 1:size(ff_images,2)
        temp_y = squeeze(ff_images(i,j,:));
        temp_y = double(temp_y);
        ind = find(itime<4);
        A = ones(numel(itime(ind)),2);
        A(:,1) = itime(ind);
        B = temp_y(ind);
        [x,stdx,mse] = lscov(A,B);
        gain(i,j) = x(1);
        offset(i,j) = x(2);
        gain_std(i,j) = stdx(1);
        offset_std(i,j) = stdx(2);
        chisq(i,j) = mse;
        %fit vs median response
    end
    %print progress
    if mod(i,10) ==0
        fprintf(1,backstr);
        str = ['Finished with row ' num2str(i) ' Of ' num2str(size(ff_images,1))];
        fprintf(1,str);
        backstr = [];
        while numel(backstr) < 2*numel(str)
            backstr = [backstr '\b'];
        end
    end            

end
