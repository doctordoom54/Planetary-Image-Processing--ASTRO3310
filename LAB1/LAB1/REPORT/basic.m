clc
clf
close all
clear
img = fitsread('imp.fits');
headerinfo = fitsheader('imp.fits');
minval = min(img);
size(minval);
size(img(:));
maxval = max(img(:));
%make the figure and place its handle in the variable h_fig1
h_fig1 = figure();
movegui(h_fig1,'center')
fig_pos = get(h_fig1,'Position');
% return the figures position on the screen
fig_pos(3:4) = [612,768];%set the figure size in display pixels to be the size (or at lest aspect ratio) of our image
set(h_fig1,'Position',fig_pos,'Name','Martian Surface');%sets the figure size/position and puts "Martian Surface" in the title bar
h_img1 = imagesc(img);%displays the image and returns the handle for the image object
h_ax1 = gca;%gets the handle for the current axis
colormap gray;%set the colormap to grayscale (default is a blue to red jet colormap)
h_cbar1 = colorbar;%put up a colorbar to show the mapping between grayscale and pixel value
axis equal;%automatically scales the axis to be equant
set(h_ax1,'XLim',[1,612],'YLim',[1,768]);
set([h_ax1,h_cbar1],'FontSize',14,'FontWeight','Bold');%makes the tick marks and tick labels easier to read
h_title1 = title('Martian Surface (MPF)','FontSize',18,'FontWeight','Bold');%gives the axes a Title
set(h_ax1,'Position',[0.05,0.05,0.9,0.9]);%expand the size of the axis in the figure to remove blank space
set(h_ax1,'YDir','Normal'); %change Y-Direction from "Reverse" to "Normal".
%%
%Stretching
if ishandle(h_fig1)
h_fig2 = figure('Position',get(h_fig1,'Position')); %make a new figure the same size and location as the first
else
h_fig2 = figure(); % make a new figure
end
h_img2 = imagesc(img,[0,0.15]); %display image stretched from 0 to 0.15
colormap gray; %set the colormap
h_cbar2 = colorbar; %open the colorbar
h_ax2 = gca; %save the handle for the new axis
set(h_ax2,'YDir','Normal'); %flip the image as before
axis equal; %automatically scales the axis to be equant
set(h_ax2,'XLim',[1,612],'YLim',[1,768]);
set([h_ax2,h_cbar2],'FontSize',14,'FontWeight','Bold'); %makes the tick marks and tick labels easier to read
h_title2 = title('Martian Surface (MPF Stretched Image)','FontSize',18,'FontWeight','Bold'); %add a title
if ishandle(h_ax1)%
set(h_ax2,'Position',get(h_ax1,'Position'));%if you still have the other figure open, set the axis to the same size
end
%%
[x,y] = ginput(4);%grabs four positions on the image by 
% left-clicking on them (if you did not include a #, 
% the function would continue grabbing values until you right-clicked your mouse).
%%
ind = sub2ind(size(img),round(y),round(x));%turns the line/sample coordinates into
% 1D band sequential indices for addressing.
z = img(ind);%address the image using 1D band sequential indices




