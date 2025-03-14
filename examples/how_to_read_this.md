fort 35 also contains some fancy boundary condition stuff, might be useful for later


deze case is van alcator c-mod, de rand klopt dus effectief


strata zijn bronnen!


fort 35: index neighborlist surface model... laatste vier zijn de koppeling met de plasma
1. cell id in plasma (-1 indien geen plasma achtergrond)
2. mapping van face 1 in triangle mesh naar nummer van de face in de plasma mesh
3. mapping van face 2 in triangle mesh naar nummer van de face in de plasma mesh
4. mapping van face 3 in triangle mesh naar nummer van de face in de plasma mesh
-> indien negatieve face index

het kan zijn dat 1 plasma face verdeeld wordt over meerede triangle faces


b2fgmtry:
    nCi -> aantal interne cellen
    nCg -> aantal guard cellen
    nCv -> totaal aantal cellen
    ...

    cvFcP -> eerste twee letters: waar, derde en vierde: wat, P pointer
    cvFC -> de faces van alle cellen
    Fccv -> de cellen voor alle faces
    fcvx -> voor elke face de vertices
    cxVx -> voor elke cell de vertices
    ft... -> flux tubes, niet relevant
    fclbl -> face labels, niet per se nodig, eventueel voor complexe surface models
    cvBb -> magnetisch veld per cell (poloidaal, radial, toroidaal en totaal, dus 4 per cell)
    cvEb -> unit vector van magnetisch veld in carthesische coordinaten
    cvX -> cell center x coord
    cvY -> cell center y coord
    cvZ -> ...
    cvSz -> oppervlakte van cel
    cvHz -> oppervlakte van loodrechte deel
    cvVol -> volume (ook gelijk aan cvSz * cvHz)
    fcs -> oppervlakte van face
    qalf, qbeta -> hoeken tussen etheta en normaal van face



b2fstati:
    ncv -> aantal cellen in plasma domein (inclusief guard cells, ik heb normaal enkel interne nodig)
    nfc -> aantal faces
    ns -> aantal species (in dit geval zijn dat atomen en ionen, ik heb enkel ionen nodig)

    zamin, zamax -> charge
    zn -> nuclear charge
    am -> atoom massa
    na -> dichtheid, op volgorde van de cellen en dan species
    ne -> electron dichtheid (niet op species)
    ua -> parallele snelheid
    uadia -> niet nodig
    te -> electron temperatuur (voor beide species)
    ti -> ion temperatuur ()
    tn -> neutrale temperatuur (hier niet zo relevant)
    po -> potential (niet nodig, eventueel voor randvoorwaarde)
    fna -> fluxen (flux voor elke face, zowel poloidaal als radiaal, kan belangrijk zijn voor de source -> wall recombination wordt bepaald door sheath model (zie eirene docs) )
    al de rest zou niet nodig moeten zijn

