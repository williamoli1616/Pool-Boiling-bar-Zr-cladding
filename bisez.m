function [xvect,xdif,fx,it]=bisez(a,b,nmax,toll,fun)
%
% [xvect,xdif,fx,it]=bisez(a,b,nmax,toll,fun) 
%
% Metodo di bisezione per la risoluzione
% dell'equazione non lineare f(x)=0
%
% Parametri di ingresso:
%
% a,b       Estremi intervallo di ricerca radice
% nmax      Numero massimo di iterazioni
% toll      Tolleranza sul test d'arresto
% fun       Function handle contenente la funzione
%
% Parametri di uscita:
%
% xvect     Vett. contenente tutte le iterate
%           calcolate (l'ultima componente e' la soluzione)
% xdif      Vett. contenente l'avanzamento tra due iterate
% fx        Vett. contenente le valutazioni di 'fun' in 'xvect'
% it        Iterazioni effettuate


% parametri di uscita
it=-1; % il primo passo sara' it=0 
xvect=[];
fx=[]; 
xdif=[];

% finche' non comincio, non so quantificare l'errore
% quindi impongo un valore fittizio all'errore
% che mi obblighi ad eseguire almeno un passo di bisezione
err=toll+1; 


% Se l'errore e' superiore alla tolleranza
% e non e' stato raggiunto il numero massimo di passi
% allora devo procedere con la bisezione
while (it < nmax && err > toll)

    %stima dello zero
    x = (b+a)/2;
    % errore
    if (fun(x) == 0)
       err=0;
    else
       err=abs(b-a)/2; 
    end
    
    % aggiornamento parametri di uscita
    xdif=[xdif;err];
    xvect=[xvect;x]; 
    fx=[fx;fun(x)];
    it=it+1;

    % scelta del nuovo estremo per l'eventuale ciclo successivo
    if (fun(x)*fun(a) < 0), 
          b=x; 
    else 
          a=x; 
    end;      
    
end; 

% Stampa sullo schermo
if (it<nmax)
    fprintf(' Convergenza al passo k : %d \n',it);
else
    fprintf(' E` stato raggiunto il numero massimo di passi k : %d \n',it);
end
fprintf(' Radice calcolata     : %-12.8f \n',xvect(it+1));

return

