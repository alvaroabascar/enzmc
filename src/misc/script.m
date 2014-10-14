load ana
load num

i = figure(1);

set(i, 'paperunits', 'inches')
set(i, 'paperorientation', 'portrait')
set(i, 'papersize', [3, 4])
set(i, 'paperposition', [0,0,4,3])

errors = (num - ana) ./ ana * 100;
plot([1:length(ana)], errors, '-or', 'linewidth', 3)

legend(['analytical', 'numerical'])

xlabel 'Km'
ylabel 'Relative error %'

print('-dpng', 'ana_num.png')
