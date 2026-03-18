#######################################################################################################
#   tSPL rule
NewRulesFor(MDDFT, rec(
    MDDFT_tSPL_RowCol := rec(
        info := "tSPL MDDFT_n -> MDDFT_n/d, MDDFT_d",

        applicable := (self, t) >> Length(t.params[1]) > 1,
        freedoms := t -> [ [1..Length(t.params[1])-1] ],

        child := (t, fr) -> let(
            newdims := SplitAt(t.params[1], fr[1]),
            rot := t.params[2],
            [ TTensor(
                MDDFT(newdims[1], rot),
                MDDFT(newdims[2], rot)
            ).withTags(t.getTags())]
        ),

        apply := (t, C, Nonterms) -> C[1],
        switch := false
    ),
    
    MDDFT_tSPL_Pease := rec(
        info := "tSPL MDDFT_n -> (I x DFT)L",

        applicable := (self, t) >> Length(t.params[1]) > 1,

        children := (self, nt) >> let(
            Ns := nt.params[1],
            N := Product(Ns),
            k := nt.params[2],
            [[ TCompose(List(Reversed(Ns), n->TTensorI(DFT(n, k), N/n, AVec, APar))).withTags(nt.getTags()) ]]
        ),

        apply := (t, C, Nonterms) -> C[1],
        switch := false
    )
));
