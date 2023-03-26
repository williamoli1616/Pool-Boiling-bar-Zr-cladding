function [xvect,it] = newton(x0,nmax,toll,fun,dfun)

% [xvect,it] = newton(x0,nmax,toll,fun,dfun) 
%
% Metodo di Newton per la ricerca degli zeri della
% funzione fun. Test d'arresto basato sul controllo
% della differenza tra due iterate successive.
%
% Parametri di ingresso:
%
% x0         Punto di partenza
% nmax       Numero massimo di iterazioni
% toll       Tolleranza sul test d'arresto
% fun dfun   Function handle contenenti la funzione e la sua derivata
%
% Parametri di uscita:
%
% xvect      Vett. contenente tutte le iterate calcolate
%            (l'ultima componente e' la soluzione)
% it         Iterazioni effettuate

err = toll+1;
it = 0;
xvect = [x0];

while (it< nmax && err>= toll)
    xv = xvect(end);
   if dfun(xv) == 0
      disp(' Arresto per azzeramento di dfun');
      it = it + 1;
      break
   else
      xn = xv - fun(xv) / dfun(xv);
      err = abs(xn - xv);
      xvect = [xvect; xn];
      it = it + 1;
   end
end

fprintf(' \n Numero di Iterazioni : %d \n',it);
fprintf(' Zero calcolato       : %-12.13f \n',xvect(end));

return

