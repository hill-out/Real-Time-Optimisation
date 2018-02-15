function [dF] = massBalance(c,Fin,k)

Fout = sum(Fin);

dF = Fin + c*Fout + k*c(1)*c(2);

end