needsPackage "InvariantRing";

n = 3; -- Anzahl der Knoten
sVars = sub(n*(n+1)/2, ZZ); -- Anzahl der Einträge der Matrix Sigma

load "colGraphs.m2";

R = QQ[for i in 1..(sVars-n) list l_i];
tries = 20; -- Anzahl der zufälligen Sigma-Matrizen, die zur Überprüfung der Identifizierbarkeit verwendet werden

-- groups() generiert alle Partitionen mit k Elementen

groups = k -> (
    if k == 1 then return {{0}};
    partsOld = groups(k-1);
    
    partsNew = {};
    for p in partsOld do (
        for c in 0..(k-2) do if p_(c) == c then partsNew = append(partsNew, append(p, c));
        partsNew = append(partsNew, append(p, k-1));
    );
    return partsNew;
);

gs = groups(sVars-n+1); -- alle Partitionen der rechten oberen Hälfte von Omega, die Kanten der 0-Partition sind nicht im Graph enthalten

-- Übersetzt die Partitionen in Kantenfärbungen. Die Eingabe p={1,0,0,1,5,0,0,8,0,10} wird beispielsweise in die Matrix  
-- (-  1  0  0  1)
-- (1  -  5  0  0)
-- (0  5  -  8  0)
-- (0  0  8  - 10)
-- (1  0  0  10 -)
-- umgewandelt. Die Farbe einer Kante i <-> j ist dann über den Eintrag an der Stelle (i,j) festgelegt. Der Graph enthält eine Kante genau dann, wenn der 
-- entsprechende Eintrag ungleich Null ist.

gs = for g in gs list (
    vals = mutableMatrix map(ZZ^n,n,(i,j)->0);
    ind = 1;
    for i in 0..(n-1) do for j in (i+1)..(n-1) do (
        if g_ind != 0 then vals_(i,j) = g_ind;
        ind = ind + 1;
    );
    vals = matrix (transpose(vals) + vals)
);

-- generiert alle möglichen (gerichtete Struktur, Kantenfärbung)-Tupel

perms = for perm in permutations n list permutationMatrix concatenate(apply(perm,x -> toString (x+1)));
configs = new MutableHashTable from for L in Ls list L => gs;
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

-- result_"special" enthält diejenigen Graphenkonfigurationen, deren Identifizierbarkeit durch die Kantenfärbung verhindert wird. Für n=4 ist diese Liste leer.

result = iterateReps(reps);