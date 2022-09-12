needsPackage "InvariantRing";

n = 5; -- Anzahl der Knoten
sVars = sub(n*(n+1)/2, ZZ); -- Anzahl der Einträge der Matrix Sigma 

load "colGraphs.m2";

-- Eingabe der betrachteten Kantenfärbung p. Die Eingabe p={1,0,0,1,5,0,0,8,0,10} wird beispielsweise in die Matrix  
-- (-  1  0  0  1)
-- (1  -  5  0  0)
-- (0  5  -  8  0)
-- (0  0  8  - 10)
-- (1  0  0  10 -)
-- umgewandelt. Die Farbe einer Kante i <-> j ist dann über den Eintrag an der Stelle (i,j) festgelegt. Der Graph enthält eine Kante genau dann, wenn der 
-- entsprechende Eintrag ungleich Null ist.

p = {1,0,0,1,5,0,0,8,0,10};
vals = mutableMatrix map(ZZ^n,n,(i,j)->0);
ind = 0;
for i in 0..(n-1) do for j in (i+1)..(n-1) do (
    if p_ind != 0 then vals_(i,j) = p_ind;
    ind = ind + 1;
);
vals = matrix (transpose(vals) + vals);

tries = 20; -- Anzahl der zufälligen Sigma-Matrizen, die zur Überprüfung der Identifizierbarkeit verwendet werden
R = QQ[for i in 1..(sVars-n) list l_i];

-- generiert alle möglichen (gerichtete Struktur, Kantenfärbung)-Tupel

perms = for perm in permutations n list permutationMatrix concatenate(apply(perm,x -> toString (x+1)));
permVals = unique for perm in perms list toStandardForm(transpose(perm)*vals*perm);
configs = new MutableHashTable from for L in Ls list L => permVals;
reps = {};

count = 0;

-- selektiert Repräsentanten und verwirft Tupel äquivalenter Graphen

for L in Ls do (
    while length configs#L > 0 do (
        count = count + 1;
        print(count);
        rep = (L,configs#L_0);
        reps = append(reps,rep);
        configs#L = delete(rep#1,configs#L);

        for perm in perms do (
            LNew = transpose(perm)*(rep#0)*perm;
            if configs#?LNew then configs#LNew = delete(toStandardForm(transpose(perm)*(rep#1)*perm),configs#LNew);
        );
    );
);

result = iterateReps(reps);

-- generiert eine Ausgabe der Anzahl identifizierbarer Graphen für jede Zahl gerichteter Kanten

f = "compAdjacentNew.csv" << "nDir, total, col, noCol" << endl;
for i in 0..(sVars-n) do (
    f << i << ", " << result#"total"#i << ", " << result#"col"#i << ", " << result#"noCol"#i << endl;
);
f << close;