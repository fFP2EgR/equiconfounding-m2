-- checkGb() überprüft das durch die Polynome in b festgelegte Modell auf Identifizierbarkeit
-- Inputs: 
--  b - erzeugende Polynome des Verschwindungsideals
--  numDir - Anzahl gerichteter Parameter
-- Output: true falls das Modell rational identifizierbar ist, sonst false

checkGb = (b,numDir) -> (
    checkList = for i in 1..numDir list false;
    checkList = new MutableList from checkList;
    
    for e in b do (
        if (degree e)_0 != 1 then continue;

        lm = leadTerm e;
        for i in 1..numDir do if lm % ideal(l_i) == 0 then (
            checkList#(i-1) = true;
            break;
        );
    );
    return all(checkList,x -> x);
);

-- toStandardForm() erhält eine Matrix O, welche eine Kantenfärbung kodiert und nummeriert diese um (das ermöglicht später den Vergleich mit anderen Kantenfärbungen).
-- Input:
--  O - Matrix mit den Einträgen c'(i,j), wobei die Farben eine beliebige Reihenfolge haben können
-- Output: Matrix mit Einträgen c(i,j), angeordnet in strikt aufsteigender Reihenfolge 

toStandardForm = O -> (
    subs = new MutableList from for i in 0..(sVars-n) list 0;
    ind = 0;
    for i in 0..(n-1) do for j in (i+1)..(n-1) do (
        ind = ind + 1;
        if O_(i,j) == 0 or subs#(O_(i,j)) != 0 then continue;
        subs#(O_(i,j)) = ind;
    );

    ONew = mutableMatrix O;
    for i in 0..(n-1) do for j in 0..(n-1) do ONew_(i,j) = subs#(O_(i,j));

    return matrix ONew;
);

-- genSigma() erzeugt eine zufällige Matrix Sigma, welche die durch p und L gegebenen Einschränkungen erfüllt.
-- Input:
--  p - kodierte Kantenfärbung
--  L - Matrix der Form (I - Lambda) 
-- Output: zufälliges Sigma

genSigma = (p,L) -> (
    subArgs = for i in 1..(sVars-n) list l_i => random(-50,50)/10;
    LSub = sub(L, subArgs);
    
    O = mutableMatrix map(QQ^n,n,(i,j)->if i == j then (random(50)+1)/20 else 0);
    OVals = {};
    for i in 1..(sVars-n) do (
        if p_i == 0 then OVals = append(OVals, 0/1)
        else if p_i == i then OVals = append(OVals, random(-50,50)/10)
        else OVals = append(OVals, OVals_(p_i-1));
    );
    
    count = 0;
    for i in 0..(n-1) do for j in (i+1)..(n-1) do (
        O_(i,j) = OVals_count;
        count = count + 1;
    );
    O = matrix O;
    O = (transpose O) + O;
    
    return inverse(transpose(LSub))*O*inverse(LSub);
);

-- checkId() generiert die erzeugenden Polynome des Verschwindungsideals und überprüft Identifizierbarkeit mit zufälligen Sigma in mehreren Durchläufen
-- Input:
--  L - Matrix der Form (I - Lambda) 
--  p - kodierte Kantenfärbung
--  numDir - Anzahl gerichteter Parameter
-- Output: true falls das Modell für ausreichend viele (>tries/2) Sigma rational identifizierbar ist, sonst false 

checkId = (L,p,numDir) -> (
    rels = mutableMatrix id_(R^(sVars-n));
    for i in 1..(sVars-n) do if p_i > 0 then rels_(i-1, p_i-1) = rels_(i-1, p_i-1) - 1; -- rels enthält lineare Abhängigkeiten zwischen den Einträgen von O
    rels = matrix rels;

    numId = 0;

    for t in 1..tries do (
        S = genSigma(p,L);
        O = transpose(L)*S*L;
            
        vec = {};
        for i in 0..(n-1) do for j in (i+1)..(n-1) do vec = append(vec, O_(i,j));
        vec = transpose matrix {vec};
        b = gens gb ideal(rels*vec);

        if checkGb((entries b)_0,numDir) then numId = numId + 1;
    );

    return numId > tries / 2;
);

-- Ls enthält Matrizen für sämtliche gerichtete Strukturen über n beschrifteten Variablen.
-- Stellvertretend für die Paramteter enthalten die Ls Einsen an den passenden Stellen

Ls = for i in 0..(2^(sVars-n)-1) list (
    L = mutableMatrix id_(ZZ^n);
    ind = 0;
    for u in 0..(n-1) do for v in (u+1)..(n-1) do (
        L_(u,v) = floor(i/2^ind) % 2;
        ind = ind + 1;
    );
    matrix L
);

-- iterateReps() erhält eine Liste mit Repräsentanten der Graphen-Äquivalenzklassen und überprüft diese auf Identifizierbarkeit, jeweils mit und ohne Kantenfärbung
-- Input:
--  reps - Liste mit Repräsentanten
-- Output: HashTable mit den unten beschriebenen Informationen

iterateReps = reps -> (
    col = new MutableList from for i in 0..(sVars-n) list 0; -- enthält für jede Zahl gerichteter Kanten die Anzahl identifizierbarer Graphen mit Kantenfärbung
    noCol = new MutableList from for i in 0..(sVars-n) list 0; -- enthält für jede Zahl gerichteter Kanten die Anzahl identifizierbarer Graphen ohne Kantenfärbung
    total = new MutableList from for i in 0..(sVars-n) list 0; -- enthält für jede Zahl gerichteter Kanten die Anzahl unterschiedlicher Graphen 

    special = {}; -- Graphen, deren Identifizierbarkeit durch eine Kantenfärbung eineschränkt wird
    len = length reps;
    ind = 1;

    for rep in reps do (
        print(toString ind | "/" | toString len);
        ind = ind + 1;

        numDir = 0;
        for i in 0..(n-1) do for j in (i+1)..(n-1) do numDir = numDir + rep#0_(i,j);
        total#numDir = total#numDir + 1;
        L = mutableMatrix id_(R^n);
        count = 1;
        for i in 0..(n-1) do for j in (i+1)..(n-1) do (
            if rep#0_(i,j) == 1 then (
                L_(i,j) = -l_count;
                count = count + 1;
            ) else (
                L_(i,j) = 0;
            );
        );
        L = matrix L;

        noColId = false;

        p = {0};
        count = 1;
        for i in 0..(n-1) do for j in (i+1)..(n-1) do (                         -- kodiert eine Kantenfärbung mit ausschließlich atomaren Komponenten, d.h. den ungefärbten Graphen
            if rep#1_(i,j) != 0 then p = append(p,count) else p = append(p,0);
            count = count + 1;
        );
        if checkId(L,p,numDir) then (
            noCol#numDir = noCol#numDir + 1;
            noColId = true;
        );

        p = {0};
        for i in 0..(n-1) do for j in (i+1)..(n-1) do if rep#1_(i,j) == 0 then p = append(p,0) else p = append(p,rep#1_(i,j)); -- kodiert die Kantenfärbung in p
        if checkId(L,p,numDir) then (
            col#numDir = col#numDir + 1;
        ) else if noColId then special = append(special,(L,p));
    );

    return new HashTable from {"col" => col, "noCol" => noCol, "total" => total, "special" => special};
);