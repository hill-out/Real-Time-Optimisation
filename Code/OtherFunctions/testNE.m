%NEtest

u0 = [3.8446,10.21,103.2121];
dPara = [1,1,1];
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];

dO.du = [1.2429,2,8.7638];

b0 = openModelPara(u0, dPara, xGuess);
b = openModelPara(u0+dO.du, dPara, xGuess);

a0 = finDiff(@(u)phiFun(u,openModelPara(u, dPara, xGuess)), u0, 0.00001)';
a = finDiff(@(u)phiFun(u,openModelPara(u, dPara, xGuess)), u0+dO.du, 0.00001)';
da = a-a0;

dO.dC = b-b0;

da2 = NEgradPara(u0+dO.du/2,dPara,dO);

e = 100*(a-(da2.dphidu'+a0))./a