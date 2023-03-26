function [t_h,u_h]=eulero_avanti(f,t0,t_max,y_0,h)

% [t_h,u_h]=eulero_avanti(f,t0,t_max,y_0,h)
%
% Risolve il problema di Cauchy 
%
% y'=f(t,y)
% y(0)=y_0
%
% utilizzando il metodo di Eulero in Avanti 
% (Eulero Esplicito) : u^(n+1)=u^n+h*f^n
%
% L'equazione differenziale ordinaria può essere in generale vettoriale 
% (y e f(t,y) appartengono a R^d, dove d=1,2,...), se d = 1 si ottiene il
% caso scalare.
%
% Input:
% -> f: function che descrive il problema di Cauchy (dichiarata tramite @). 
%       Deve ricevere in ingresso due argomenti: f=f(t,y), di cui y è in
%       generale vettore di lunghezza d.
% -> t0, t_max: estremi dell' intervallo temporale di soluzione 
% -> y_0: il dato iniziale del problema di Cauchy (un vettore colonna di
%       lunghezza d)
% -> h: l'ampiezza del passo di discretizzazione temporale.
% ATTENZIONE: si controlli che l'output di f e il dato y_0 siano entrambi
% vettori colonna della stessa lunghezza!
%
% Output:
% -> t_h: vettore degli istanti in cui si calcola la soluzione discreta
% -> u_h: la soluzione discreta calcolata nei nodi temporali t_h (in ogni 
%       istante, la soluzione è un vettore colonna; u_h è perciò una matrice
%       rettangolare di dimensioni d x N_istanti)

% vettore degli istanti in cui risolvo la eq diff ord
t_h=t0:h:t_max;

% inizializzo il vettore che conterra' la soluzione discreta
N_istanti=length(t_h);
d = length(y_0);
u_h=zeros(d,N_istanti);

% ciclo iterativo che calcola u_(n+1)=u_n+h*f_n
u_h(:,1)=y_0;

for it=2:N_istanti
    u_old=u_h(:,it-1);
    u_h(:,it)=u_old+h*f(t_h(it-1),u_old);
end

