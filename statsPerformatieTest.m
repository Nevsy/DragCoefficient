% metingen test 1 en 2
metingen1 = [7, 6, 10, 7, 3, 3, 3, 3, 3, 0];
metingen2 = [3, 3, 3, 3, 0, 0, 0, 0, 0, 0];

std1 = std(metingen1);
std2 = std(metingen2);
mu1 = sum(metingen1)/length(metingen1);
mu2 = sum(metingen2)/length(metingen2);