
# Copyright (c) 2018-2021, Carnegie Mellon University
# See LICENSE for details


Class(fma, Exp, rec(
    computeType := self >> UnifyTypesL(self.args),
    toMul := self >> self.args[1] + self.args[2] * self.args[3],
    ev := self >> self.args[1].ev() + self.args[2].ev() * self.args[3].ev()));

Class(fms, Exp, rec(
    computeType := self >> UnifyTypesL(self.args),
    toMul := self >> self.args[1] - self.args[2] * self.args[3],
    ev := self >> self.args[1].ev() - self.args[2].ev() * self.args[3].ev()));

Class(nfma, Exp, rec(
    computeType := self >> UnifyTypesL(self.args),
    toMul := self >> -self.args[1] -+self.args[2] * self.args[3],
    ev := self >> -self.args[1].ev() + self.args[2].ev() * self.args[3].ev()));

RewriteRules(RulesStrengthReduce, rec(
   fma_Xab := Rule([@(0,[fma,fms,nfma]), _vtrivial, @, @], e -> e.toMul()), 
   fma_aXb := Rule([@(0,[fma,fms,nfma]), @, _vtrivial, @], e -> e.toMul()), 
   fma_abX := Rule([@(0,[fma,fms,nfma]), @, @, _vtrivial], e -> e.toMul()), 
));

notMul := e -> ObjId(e)<>mul;

#NOTE - shouldn't have to handle arrays here
_IsAVector  := (var) >> (IsVecT(var.t) or IsArray(var.t));

_isNonZero  := i -> i.ev() <> 0;

_isDivOk    := e -> When(_IsAVector(e), ForAll(e.v, i -> _isNonZero(i)), _isNonZero(e));

_isNotEqual := function(a, b)
   # Why did we check for vectors?
   #if _IsAVector(a) or _IsAVector(b) then
   #   return true;
   #fi;
   return (a<>b);
end;

_divide := function(numerator, denominator)
# Both are scalars
   local type, vLength, nVal, dVal;
   if not (_IsAVector(numerator) or _IsAVector(denominator)) then
      return numerator/denominator;
   fi;

   if _IsAVector(numerator) and not _IsAVector(denominator) then
      type     := numerator.t;
      vLength  := numerator.t.size;
      nVal     := numerator.v;
      return(Value(type, List([1..vLength], i->nVal[i]/denominator.v)));
   fi;

   if not _IsAVector(numerator) and _IsAVector(denominator) then
      type     := denominator.t;
      vLength  := denominator.t.size;
      dVal     := denominator.v;
      return(Value(type, List([1..vLength], i->numerator.v/dVal[i])));
   fi;

# Both are vectors. Assume lengths are same.
   type     := numerator.t;
   vLength  := numerator.t.size;
   nVal     := numerator.v;
   dVal     := denominator.v;
   return(Value(type, List([1..vLength], i->nVal[i]/dVal[i])));
end;

_def := function(list, t, exp) 
    local cmd;
    cmd := assign(var.fresh_t("s",t), exp);
    cmd.loc.def := cmd;
    Add(list, cmd);
    return cmd.loc;
end;

FMA := function(code)
    local c, cmds, dd, exp, res, supp, ch;
    res := [];

    cmds := code.cmds;

    for c in cmds do
        # NOTE: needs to check for our FP datatype, not 'uncheck' for TPtr, TInt etc.
        # NOTE: determine how pointer arithmetic works with FMAs on various archs.
        #       For now, avoid FMA'ing index computations. 
        if IsBound(c.loc) and not IsPtrT(c.loc.t) and not IsIntT(c.loc.t) then
            dd := DefCollect(c,2); # definition, depth=2

            # ch is like children of DefCollect(c,1)
            if IsBound(c.exp.args) then ch := c.exp.args; else ch := []; fi; 
            supp := [];

            exp := MatchSubst(dd, [
                [ [mul, @(1,Value), [mul, @(2,Value), @(3)]], 
                    e -> mul(@(1).val*@(2).val, @(3).val), "mul mul" ],
 
                [ [add, [mul, @(1, Value), @(2)],
                        [mul, @(3, Value, e->_isDivOk(e) or _isDivOk(@1.val)), @(4)]], 
                #NOTE: Determine small/large if scalar (we now only determine
                #       zero/nonzero if vector)
                    e -> let(
                        sm := When(_isDivOk(@(1).val), @(1).val, @(3).val), 
                        lg := When(_isDivOk(@(1).val), @(3).val, @(1).val),
                        Cond(_isNotEqual(@(1).val, @(3).val),
                            Cond(sm.v=@(1).val.v,
                                # propagate (divide by) the smaller value
                                sm * _def(supp, c.loc.t, fma(@(2).val, 
                                    _divide(@(3).val, @(1).val), @(4).val)), 
                                sm * _def(supp, c.loc.t, fma(@(4).val, 
                                    _divide(@(1).val, @(3).val), @(2).val))
                            ), 
                            sm * _def(supp, c.loc.t, add(@(2).val, @(4).val))
                        )),
                    "fma:ab+cd"],
                
                [ [sub, [mul, @(1, Value), @(2)], 
                        [mul, @(3, Value, e->_isDivOk(e) or _isDivOk(@1.val)), @(4)]], 
                # NOTE: Determine small/large if scalar (we now only determine
                #        zero/nonzero if vector)
                    e -> let(
                        sm := When(_isDivOk(@(1).val), @(1).val, @(3).val),
                        lg := When(_isDivOk(@(1).val), @(3).val, @(1).val),
                        Cond(_isNotEqual(@(1).val, @(3).val),
                            Cond(sm.v=@(1).val.v,
                                sm * _def(supp, c.loc.t, fms(@(2).val, 
                                    _divide(@(3).val, @(1).val), @(4).val)), 
                                sm * _def(supp, c.loc.t, nfma(@(4).val, 
                                    _divide(@(1).val, @(3).val), @(2).val))
                            ), 
                            sm * _def(supp, c.loc.t, sub(@(2).val, @(4).val))
                        )),
                    "fma:ab-cd"],

                [ [add, [mul, @(1), @(2)], @(3)], 
                    e -> fma(ch[2], @(1).val, @(2).val), "fma"],
                [ [add, @(3), [mul, @(1), @(2)]], 
                    e -> fma(ch[1], @(1).val, @(2).val), "fma2"],

                [ [sub, [mul, @(1), @(2)], @(3)], 
                    e -> nfma(ch[2], @(1).val, @(2).val), "nfma" ],
                [ [sub, @(3), [mul, @(1), @(2)]], 
                    e -> fms(ch[1], @(1).val, @(2).val), "fms" ],

            ]);

            if not Same(dd, exp) then
                #Print("-----------------------------------------------------\n");
                #Print("c:    ", c  ,  "\n-----\n"); Print("dd:   ", dd ,  "\n-----\n");
                #Print("exp:  ", exp,  "\n-----\n"); Print("supp: ", supp, "\n-----\n");
                #Print("-----------------------------------------------------\n");
                c.exp := exp; 
                c.loc.pred := ArgsExp(exp);
                Append(res, supp);
                Add(res, c);
            else 
                Add(res, c);
            fi;
        else
            Add(res, c);
        fi;
    od;

    return chain(res);
end;

DoFMA := function(code)
   code := BinSplit(code);
   MarkDefUse(code);
   code := FMA(code);
   MarkDefUse(code);
   #return CopyPropagate(code); # Do this in the compileStrategy
   return(code);
end;

