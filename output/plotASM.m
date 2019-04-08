close all; clear all; clc; warning off
%% Plotting of output data from ASM
A=importdata('meandata.dat');
data = A.data;
x = data(:,1);
y = data(:,2);
ugrid = data(:,3);
vgrid = data(:,4);

quiver(x,y,ugrid,vgrid)
axis image


%% Make gif of bubbly flow over time
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
iterations = 399;
data_points = num2str([0:iterations].','0%04d')
formatSpec = "d%s.dat";

for n = 1:length(data_points(:,1))
    data_file = sprintf(formatSpec,data_points(n,:));
    A=importdata(data_file);
    data = A.data;
    x = data(:,1);
    y = data(:,2);
    ugrid = data(:,3);
    vgrid = data(:,4);
    quiver(x,y,ugrid,vgrid)
    axis image
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end
    clear data_file x y ugrid vgrid
  end

