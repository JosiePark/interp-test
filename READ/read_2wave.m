clear all;

cd '/home/josiepark/Project/PhD/DATA/INTERP/TEST_CASES/TWO_WAVE/NO_RK4';
ncid = fopen('2wave.dat','rb');

ii = 512;
jj= 512;

x = 1:ii;
y = 1:jj;

[X,Y] = meshgrid(x,y);

fread(ncid,1,'int32');
psi = fread(ncid,[ii,jj],'float64');


ncid1 = fopen('u.dat','rb');
fread(ncid1,1,'int32');
u = fread(ncid1,[ii,jj],'float64');
ncid2 = fopen('v.dat','rb');
fread(ncid2,1,'int32');
v = fread(ncid2,[ii,jj],'float64');

cd('/home/josiepark/Project/PhD/FIGURES/INTERP/TEST_CASES/2WAVE');

fig = figure('Position',[.25,.25,2250,500]);

subplot(1,3,1)
contourf(X,Y,u);
colormap('jet');
xlabel('x (grid point)');
ylabel('y (grid point)');
axis square;
colorbar
title({'u'});

subplot(1,3,2)
contourf(X,Y,v);
colormap('jet');
xlabel('x (grid point)');
ylabel('y (grid point)');
axis square;
colorbar
title({'v'});

subplot(1,3,3)
contourf(X,Y,psi);
colormap('jet');
xlabel('x (grid point)');
ylabel('y (grid point)');
axis square;
colorbar
title({'Stream function'});

print('2wave_all','-djpeg');

%%
clear all;
cd '/home/josiepark/Project/PhD/DATA/INTERP/TEST_CASES/TWO_WAVE/NO_RK4';
npoints = 1;



    
    traj_bi = ncread(sprintf('bicubic_dt_ 2160.nc'),'Trajectories');
    traj_2D = ncread(sprintf('2Dcubic_dt_ 2160.nc'),'Trajectories');
    traj_ex = ncread(sprintf('exact_dt_ 2160.nc'),'Trajectories');
    
fig = figure('Position',[.25,.25,2250,500]);
    subplot(1,3,1)
    plot(squeeze(traj_bi(1,1,:)),squeeze(traj_bi(2,1,:)));
    xlabel('x (grid point)');
    ylabel('y (grid point)');
    xlim([260,340]);
    ylim([260,340]);
    axis square;
    title(sprintf('Bicubic Interpolation'));


    subplot(1,3,2)
    plot(squeeze(traj_2D(1,1,:)),squeeze(traj_2D(2,1,:)));
    xlabel('x (grid point)');
    ylabel('y (grid point)');
    axis square;
    xlim([260,340]);
    ylim([260,340]);
    title(sprintf('2D-cubic Interpolation'));

    subplot(1,3,3)
    plot(squeeze(traj_ex(1,1,:)),squeeze(traj_ex(2,1,:)));
    xlabel('x (grid point)');
    ylabel('y (grid point)');
    axis square;
    xlim([260,340]);
    ylim([260,340]);
    title(sprintf('No Interpolation (exact)'));
    
  %% calculate how many loops were made and the deviation from initial position
eps = 0.5;
runs = 1
nt = size(traj_2D,3);
for j = 1:runs;

x0_2D(j) = traj_2D(1,1,1);
y0_2D(j) = traj_2D(2,1,1);
x0_bi(j) = traj_bi(1,1,1);
y0_bi(j) = traj_bi(2,1,1);
x0_ex(j) = traj_ex(1,1,1);
y0_ex(j) = traj_ex(2,1,1);

count = 0;
for i = 1:nt(j)-1;
    if (traj_bi(1,1,i) > x0_bi(j)) && (traj_bi(1,1,i) < x0_bi(j) + eps) && (traj_bi(1,1,i+1) < x0_bi(j)) && (traj_bi(2,1,i) <300)
        count = count + 1;
        crossing_bi(count) = traj_bi(2,1,i);
        cross_time_bi(count) = i;
    end
end


count = 0;
for i = 1:nt(j)-1;
    if (traj_2D(1,1,i) > x0_2D(j)) && (traj_2D(1,1,i) < x0_2D(j) + eps) && (traj_2D(1,1,i+1) < x0_2D(j)) && (traj_2D(2,1,i) < 300)
        count = count + 1;
        crossing_2D(count) = traj_2D(2,1,i);
        cross_time_2D(count) = i;
    end
end
count = 0;
for i = 1:nt(j)-1;
    if (traj_ex(1,1,i) > x0_ex(j)) && (traj_ex(1,1,i) < x0_ex(j) + eps) && (traj_ex(1,1,i+1) < x0_ex(j)) && (traj_ex(2,1,i) < 300)
        count = count + 1;
        crossing_ex(count) = traj_ex(2,1,i);
        cross_time_ex(count) = i;
    end
end
end

%%

for i = 1:runs;
    
    error_2D(i) = abs(crossing_2D(1) - crossing_2D(end))/crossing_2D(1);
    error_bi(i) = abs(crossing_bi(1) - crossing_bi(end))/crossing_2D(1);
    error_ex(i) = abs(crossing_ex(1) - crossing_ex(end))/crossing_ex(1);
    
end

%%

size(crossing_2D{1},2)

%%
for i = 1:runs;

loop_err_2D(i) = error_2D(i)/size(crossing_2D,2);
loop_err_ex(i) = error_ex(i)/size(crossing_ex,2);
loop_err_bi(i) = error_bi(i)/size(crossing_bi,2);

end

   
    
