%A3310 LAB2 WORKBOOK SPRING 2023
%Name: Sachin Kelkar______
%Date: xxx______

%include the Matlab commands you used to answer the LAB questions in this
%workbook. If run, the workbook should generate all figure and save all
%images/arrays relevant to the lab without any errors.

%add the LAB1 sub-directories to the Matlab path
addpath(genpath('../'));

%INSERT YOU COMMANDS BELOW (REMEMBER TO COMMENT YOUR CODE!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
f = dir("../DATA/*.IMG");
for i = 1:numel(f)
    info = pds_label_parse_v3([f(i).folder '/' f(i).name]);
    tint(i) = info.instrument_state_parms.exposureduration;
    flux(i) = info.derived_image_parms.sphereflux;
    efl(i) = info.instrument_state_parms.focallength;
    if mod(i,25) ==0
        disp(['Finished' num2str(i) 'of' num2str(numel(f))])
    end
end
%%
close all
bayer_red = imread('../DATA/bayer_mask_red.tif');
figure
imagesc(bayer_red')
ax1 = gca;
bayer_green = imread('../DATA/bayer_mask_green.tif');
figure
imagesc(bayer_green')
ax2 = gca;
linkaxes([ax1,ax2])
red = csvread('MastcamZ_Throughput_InBand_R0_r_V2.csv');
blue = csvread('MastcamZ_Throughput_InBand_R0_b_V2.csv');
green = csvread('MastcamZ_Throughput_InBand_R0_g_V2.csv');
figure 
plot(red(:,1),red(:,2),'r')
hold on 
plot(green(:,1),green(:,2),'g')
hold on
plot(blue(:,1),blue(:,2),'b')
