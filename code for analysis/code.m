%% imread
I=imread('C:\project/6X AF-4_1.jpg');
imagesc(I);
%% Locate the ROI
center = zeros(2,1);
delta = zeros(2,1);
for iDir = 1:2
    %% Project data
    proj = sum(I(:,:,1),iDir);
    %% Find edge points
    thresh = max(proj)*.5;
    i1 = find(proj >= thresh,1,'first');
    i2 = find(proj >= thresh,1,'last');
    center(iDir) = mean([i1 i2]);
    delta(iDir) = i2 - i1;
end
s = 1024;%min(delta);
ROI_X = round(center(1) + [-s/2:s/2]);
ROI_Y = round(center(2) + [-s/2:s/2]);
%%
I_crop = I(ROI_Y,ROI_X,1);
%%
imagesc(I)
hold on
plot(ROI_X,ROI_Y(1)*ones(size(ROI_X)),'Linewidth',3)
plot(ROI_X,ROI_Y(end)*ones(size(ROI_X)),'Linewidth',3)
%%
hold off
imagesc(I_crop)
%% Calculate the radial-average power spectrum
Iu = abs(fft2c(I_crop));

%%
mask = ones(size(Iu));
% mask(round(end/2),round(end/2)) = 0;
imagesc(Iu.*mask)
%%
I_rad = radial_avg(Iu,size(Iu,1));
%%
plot(I_rad/max(I_rad))
%% Find the first minimum
[pks,locs] = findpeaks(-I_rad);
i1 = locs(1)+5;
%%
% hold on
% plot(locs,-pks,'o')
%% Find the last peak
thresh = 0.02;
i2 = find(I_rad >= max(I_rad)*thresh,1,'last');
%%
ii = i1:i2;
y = I_rad(ii);
plot(ii,y)
%% Nonlinear least-squares fit
f = @(c,x) c(4) + c(1)*exp(-x.^2/(2*c(2)^2)).*abs(sin(c(3)*x.^2));
c0 = [max(y), 15, 0.01, min(y)];
plot(ii,y,ii,f(c0,ii))
%%
cost = @(c) sum( ( f(c,ii) - y).^2 );

c_fit = fminunc(@(c) cost(c),c0);

%%
plot(ii,y,ii,f(c_fit,ii))
%%
N = floor(size(Iu,1)/2);
fx = -N:N;
[fxG,fyG] = meshgrid(fx);
fR = sqrt(fxG.^2 + fyG.^2);
mask = and(fR >= i1,fR <= i2);
imagesc(mask.*Iu)
%% 2D fit
y = Iu(mask);
cost2D = @(c) sum( ( f(c,fR(mask)) - y).^2 );
c0 = [max(y)*2, 20, 0.0053, min(y)];

ii = abs(fx) <= max(fR(mask));

y_guess = f(c0,fR(mask));
guess = zeros(size(Iu));
guess(mask) = y_guess;

% imagesc(mask(ii,ii).*Iu(ii,ii))
imagesc(mask(ii,ii).*guess(ii,ii))
%%
imagesc(mask(ii,ii).*guess(ii,ii) - mask(ii,ii).*Iu(ii,ii))
%%
c_fit_2D = fminunc(@(c) cost2D(c),c0);

y_fit = f(c_fit_2D,fR(mask));
fit2D = zeros(size(Iu));
fit2D(mask) = y_guess;

imagesc(mask(ii,ii).*fit2D(ii,ii))
%% Grid-search




