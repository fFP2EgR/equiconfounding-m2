n = 5; -- Anzahl der Knoten
sVars = sub(n*(n+1)/2, ZZ); -- Anzahl der Einträge der Matrix Sigma 

load "colGraphs.m2";

-- Eingabe der betrachteten gerichteten Struktur. Die Eingabe p={1,0,0,0,1,0,0,1,0,1} wird beispielsweise in die Matrix  
-- (-  1  0  0  0)
-- (-  -  1  0  0)
-- (-  -  -  1  0)
-- (-  -  -  -  1)
-- (-  -  -  -  -)
-- umgewandelt. Der Graph enthält eine Kante i -> j genau dann, wenn der Eintrag an der Stelle (i,j) gleich Eins ist.

d = {1,0,0,0,1,0,0,1,0,1};

-- Eingabe der betrachteten Kantenfärbung p. Die Eingabe p={1,0,0,1,5,0,0,8,0,10} wird beispielsweise in die Matrix  
-- (-  1  0  0  1)
-- (1  -  5  0  0)
-- (0  5  -  8  0)
-- (0  0  8  - 10)
-- (1  0  0  10 -)
-- umgewandelt. Die Farbe einer Kante i <-> j ist dann über den Eintrag an der Stelle (i,j) festgelegt. Der Graph enthält eine Kante genau dann, wenn der 
-- entsprechende Eintrag ungleich Null ist.

p = {0,2,3,4,0,6,7,0,9,0};

R = QQ[for i in 1..(sVars-n) list l_i];

L = mutableMatrix id_(R^n);
ind = 0;
count = 1;
for i in 0..(n-1) do for j in (i+1)..(n-1) do (
    if d_ind == 1 then (
        L_(i,j) = -l_count;
        count = count+1;
    );
    ind = ind + 1; 
);

result = checkId(matrix L,join({0},p),sum(d)); -- true falls der Graph rational identifizierbar ist, sonst false