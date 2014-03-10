indel.char.scores<-function(align) {
  scores<-c()
  score<-0
  for (i in 1:nchar(align)) {
    type<-substr(align, i, i)
    if (type=='I') {
      score<-score+1
    } else if (type=='D') {
      score<-score-1
    }
    scores<-append(scores, score)
  }
  scores
}

identify.disagree.regions<-function(exp.scores, act.scores) {
  align.len<-length(exp.scores)
  starts<-c()
  ends<-c()
  open<-FALSE
  for (i in 1:align.len) {
    if (!open) {
      if (exp.scores[[i]]!=act.scores[[i]]) {
        open<-TRUE
        starts<-append(starts, i)
      }
    } else {
      if (exp.scores[[i]]==act.scores[[i]]) {
        open<-FALSE
        ends<-append(ends, i)
      }
    }
  }
  data.frame(Start=starts, End=ends)
}

align.equal<-function(exp.align, act.align) {
  comp.boolean<-FALSE
  if (nchar(exp.align)!=nchar(act.align)) {
    comp.boolean<-FALSE
  } else {
    disagreed.regions<-identify.disagree.regions(indel.char.scores(exp.align), indel.char.scores(act.align))
    if (nrow(disagreed.regions)==0) {
      comp.boolean<-TRUE
    } else {
      for (i in 1:nrow(disagreed.regions)) {
        comp<-small.segment.align.equal(substr(exp.align, disagreed.regions$Start,
                                               disagreed.regions$End),
                                        substr(act.align, disagreed.regions$Start,
                                               disagreed.regions$End))
        if (!comp) {
          comp.boolean<-FALSE
          break
        } else {
          comp.boolean<-TRUE
        }
      }      
    }
  }
  comp.boolean
}

num.changetype<-function(align, type) {
  num<-sum(attr(gregexpr(paste(type, '+', sep=''), align)[[1]], 'match.length'))
  if (num<0) num<-0
  num
}

small.segment.align.equal<-function(exp.align, act.align, tol=0.5) {
  types<-c('C', 'M', 'I', 'D')
  exp.type.counts<-mapply(function(type) {
    num.changetype(exp.align, type)
  }, types)
  act.type.counts<-mapply(function(type) {
    num.changetype(act.align, type)
  }, types)
  comp.boolean<-FALSE
  equal.map<-(exp.type.counts==act.type.counts)
  if (equal.map[['I']] & equal.map[['D']]) {
    if (equal.map[['M']] & equal.map[['C']]) {
      comp.boolean<-TRUE
    } else {
      Mtol<-abs(exp.type.counts[['M']]-act.type.counts[['M']])/max(exp.type.counts[['M']], act.type.counts[['M']])
      Ctol<-abs(exp.type.counts[['C']]-act.type.counts[['C']])/max(exp.type.counts[['C']], act.type.counts[['C']])
      comp.boolean<-(Mtol<=tol & Ctol<=tol)
    }
  } else {
    comp.boolean<-FALSE
  }
  comp.boolean
}