function dist = dist_K(tm,xm,td,xd)
dt = tm(2) - tm(1);
index = floor(td./dt);

% calc distance between model and data
distance = (xm(index,:) - xd).^2;

% MSE
dist = mean(mean(distance,2));
end