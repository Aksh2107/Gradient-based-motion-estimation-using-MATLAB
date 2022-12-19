close all
clear

close all
clear

pic = imread('hanks.big.001.png', 'png');
pic1 = imresize(pic, 0.5);
[rows, cols, chan] = size(pic1);

pic = imread('hanks.big.002.png', 'png');
pict = imresize(pic, 0.5);
pic2 = pict(1 : rows, 1 : cols, 1 : chan);

x = 1 : cols;
y = 1 : rows;
X = ones(rows, 1) * x;
Y = (1 : rows)' * ones(1, cols);

figure(1);
image(x, y, pic1);
axis image;
hold on;
figure(2);
image(x, y, pic2);
axis image;
hold on;
shg;



% Bigger block
xx = 300;
yy = 50;
B = 128;

% % Small block in the eye
% xx = 460;
% yy = 112;
% B = 32;
% 
% %Block in the sky
% xx = 710;
% yy = 250;
% B = 128;


sy = (yy : yy + B - 1);
sx = (xx : xx + B - 1);
%B = 128;
%sy = (50 : 50 + B - 1);
%sx = (300 : 300 + B - 1);
% sy = 120 : 120 + B - 1;
% sx = 10 : 10 + B - 1;
SX = ones(B, 1) * sx;
SY = sy' * ones(1, B);

block1 = double(pic1(sy, sx, 2));
block2 = double(pic2(sy, sx, 2));

figure(3);
image(block1);
colormap(gray(256));
figure(4);
image(block2);
colormap(gray(256));
shg;

figure(5);
image(drawbox(pic1, min(sy), min(sx), length(sy), length(sx)));
axis image;
axis([min(sx) - 100, max(sx) + 100, min(sy) - 100, max(sy) + 100]);
dx = 0;
dy = 0;
block_curr = double(pic2(sy, sx, 2));
block_prev = double(pic1(sy, sx, 2));
figure(6);
image([block_curr block_prev (block_curr - block_prev + 128)]);
colormap(gray(256));
hold on;
max_iters = 20;
par_hist = zeros(max_iters, 5);
for i = 1 : max_iters,
  block_prev = interp2(X, Y, double(pic1(:, :, 2)), SX + dx, SY + dy);
  [G E update] = calcgradest(block_prev, block_curr);
  dx = dx - update(1);
  dy = dy - update(2);
  par_hist(i, :) = [dx dy update(1) update(2) mean(mean(E.^2))];
  figure(5);
  image(drawbox(pic1, round(min(sy + dy)), round(min(sx + dx)), length(sy), length(sx)));
  axis image;
  axis([min(sx) - 100, max(sx) + 100, min(sy) - 100, max(sy) + 100]);
  shg
  
  figure(6);
  image([block_curr block_prev (block_curr - block_prev + 128)]);
  axis image;
  shg
  pause
 
end;

figure(7);
plot((1 : max_iters), par_hist(:, 1), 'r-x', (1 : max_iters), par_hist(:, 2), 'b-x', 'linewidth' , 1.5);
title('Covergence of motion');
legend('DX', 'DY');
xlabel('Iteration');
ylabel('Motion');
shg;

figure(8);
plot((1 : max_iters), par_hist(:, 5), 'r-x', 'linewidth' , 1.5);
title('Covergence of Motion compensated MSE');
xlabel('Iteration');
ylabel('Block MSE');
shg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G E update] = calcgradest(block_prev, block_curr)
  [gx, gy] = gradient(block_prev);
  E = block_curr - block_prev;
  G = [gx(:) gy(:)]';
  GtG = [sum(sum(gx .* gx)) sum(sum(gx .* gy)); sum(sum(gx .* gy)) sum(sum(gy .* gy))];
  update = -inv(GtG) * G * E(:); 
end

function boxpic = drawbox(pic, row, col, N, M)
  boxpic = pic;
  [rows, cols, chan] = size(pic);
  if (chan == 3)
    x = (1 : M) + col - 1;
    y = (1 : N) + row - 1;
    boxpic(row, x, 1) = 255 * ones(1, M);
    boxpic(row + 1, x, 1) = 255 * ones(1, M);
    boxpic(row + N - 1, x, 1) = 255 * ones(1, M);
    boxpic(row + N - 2, x, 1) = 255 * ones(1, M);
    
    for n = 2 : 3,
      boxpic(row, x, n) = 0 * ones(1, M);
      boxpic(row + 1, x, n) = 0 * ones(1, M);
      boxpic(row + N - 1, x, n) = 0 * ones(1, M);
      boxpic(row + N - 2, x, n) = 0 * ones(1, M);
    end;
    
    boxpic(y, col, 1) = 255 * ones(N, 1);
    boxpic(y, col + 1, 1) = 255 * ones(N, 1);
    boxpic(y, col + M - 1, 1) = 255 * ones(N, 1);
    boxpic(y, col + M - 2, 1) = 255 * ones(N, 1);
    
    for n = 2 : 3,
      boxpic(y, col, n) = 0 * ones(N, 1);
      boxpic(y, col + 1, n) = 0 * ones(N, 1);
      boxpic(y, col + M - 1, n) = 0 * ones(N, 1);
      boxpic(y, col + M - 2, n) = 0 * ones(N, 1);
    end;
    
 
  else
    fprintf('Need 3 channel image \n');
  end

end


