close all; clear all; clc; warning off
%% Specify results file
NPJ = 150;

%% Reading data from results ASM file
A=importdata('meandata.dat');
data = A.data;
x = data(:,1);
y = data(:,2);
u = data(:,3);
v = data(:,4);
alpha = data(:,8);

%% Modify datastructure and plot mean values
xgrid = meshgrid(x, NPJ);
ygrid = meshgrid(y, NPJ);
ugrid = meshgrid(u, NPJ);
vgrid = meshgrid(v, NPJ);
alphagrid = meshgrid(alpha, NPJ);

for i = size(ugrid,1)
    for j = size(ugrid,2)
        velnorm(i,j) = norm([ugrid(i,j), vgrid(i,j)]);
    end
end

pcolor(xgrid,ygrid,alphagrid);
axis('image')
colorbar

%% Make gif of bubbly flow over time
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
iterations = 399;
data_points = num2str([0:iterations].','0%04d');
formatSpec = "d%s.dat";

for n = 1:length(data_points(:,1))
    data_file = sprintf(formatSpec,data_points(n,:));
    A=importdata(data_file);
    data = A.data;
    x = data(:,1);
    y = data(:,2);
    alpha = data(:,8);
    xgrid = meshgrid(x, NPJ);
    ygrid = meshgrid(y, NPJ);
    alphagrid = meshgrid(alpha, NPJ);
    pcolor(xgrid, ygrid, alphagrid)
    axis image
    colorbar
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