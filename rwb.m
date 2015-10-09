function colors = rwb(n)

colors = [
 1 0 0
 1 1 1
 0 0 1
];

colors = flipdim(colors,1);
colors_up = zeros(n,3);


for k = 1:3
colors_up(:,k) = interp1(colors(:,k),linspace(1,length(colors),n));
end

colors = colors_up;
