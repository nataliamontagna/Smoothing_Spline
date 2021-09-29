% [ta,ya,wa] = splav(t,y,w) calcola la media dei valori che si 
% ripetono in più nodi
%
% t = nodi con ripetizioni
% y = valori osservati nei nodi con ripetizioni
% w = pesi dei valori di t
%
% ta = nodi senza ripetizioni
% ya = media delle osservazioni corrispondenti ai nodi in ta
% wa = pesi aggiornati con le molteplicità delle osservazioni ripetute

function [ta,ya,wa] = splav(t,y,w)

if nargin==0, help splav; return; end
if nargin==2, w = ones(length(t),1); end

c = size(t,1);

[t,j] = sort(t);
y = y(j);
w = w(j);

t = t(:)'; y = y(:)'; w = w(:)';

ia = [1, find(diff(t)~=0)+1];
ta = t(ia);
Na = length(ta);

for i=1:Na-1
   r = ia(i) : ia(i+1)-1;
   wa(i) = sum(w(r));
   ya(i) = sum(w(r).*y(r))/wa(i);
end

r = ia(Na) : length(t);
wa(Na) = sum(w(r));
ya(Na) = sum(w(r).*y(r))/wa(Na);

if c>1, ta=ta'; ya=ya'; wa=wa'; end