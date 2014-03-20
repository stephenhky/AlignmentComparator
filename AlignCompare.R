indel.char.scores<-function(align) {
  scores<-c()
  score<-0
  for (i in 1:nchar(align)) {
    type<-substr(align, i, i)
    delta<-ifelse(type=='I', 1, ifelse(type=='D', -1, 0))
    score<-score+delta
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
        if (i==align.len) { ends<-append(ends, i) }
      }
    } else {
      if ((exp.scores[[i]]==act.scores[[i]]) || (i==align.len)) {
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
        reg.start<-disagreed.regions$Start[[i]]
        reg.end<-disagreed.regions$End[[i]]
        comp<-small.segment.align.equal(substr(exp.align, reg.start, reg.end),
                                        substr(act.align, reg.start, reg.end))
        if (is.na(comp)) {comp.boolean<-FALSE}
        comp.boolean<-comp
        if (!comp) {break}
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
  exp.type.counts<-mapply(function(type) {num.changetype(exp.align, type)}, types)
  act.type.counts<-mapply(function(type) {num.changetype(act.align, type)}, types)
  comp.boolean<-FALSE
  equal.map<-(exp.type.counts==act.type.counts)
  if (equal.map[['I']] & equal.map[['D']]) {
    if (equal.map[['M']] & equal.map[['C']]) {
      comp.boolean<-TRUE
    } else {
      Mtol<-abs(exp.type.counts[['M']]-act.type.counts[['M']])/max(exp.type.counts[['M']], act.type.counts[['M']])
      Ctol<-abs(exp.type.counts[['C']]-act.type.counts[['C']])/max(exp.type.counts[['C']], act.type.counts[['C']])
      if (is.na(Mtol)) {
        comp.boolean<-(Ctol<=tol | abs(exp.type.counts[['C']]-act.type.counts[['C']])<=1)
      } else if (is.na(Ctol)) {
        comp.boolean<-(Mtol<=tol | abs(exp.type.counts[['M']]-act.type.counts[['M']])<=1)
      } else {
        comp.boolean<-((Mtol<=tol | abs(exp.type.counts[['M']]-act.type.counts[['M']])<=1) & (Ctol<=tol | abs(exp.type.counts[['C']]-act.type.counts[['C']])<=1)) 
      }
    }
  } else {
    comp.boolean<-FALSE
  }
  comp.boolean
}

smart.align.equal<-function(exp.align, act.align) {
  comp.boolean<-FALSE
  if (nchar(exp.align)==nchar(act.align)) {
    comp.boolean<-align.equal(exp.align, act.align)
  } else {
    longer.align<-ifelse(nchar(exp.align)>nchar(act.align), exp.align, act.align)
    shorter.align<-ifelse(nchar(exp.align)<nchar(act.align), exp.align, act.align)
    len.diff<-nchar(longer.align)-nchar(shorter.align)
    for (i in 1:len.diff) {
      if (align.equal(shorter.align, substr(longer.align, 1+i-1, nchar(longer.align)-len.diff+i-1))) {
        comp.boolean<-TRUE
        break
      }
    }
  }
  comp.boolean
}

score.align<-function(align, match.score=5, mismatch.score=-4, gap.open.score=-30, gap.extend.score=-1) {
  gap.opened<-FALSE
  score<-0
  for (i in 1:nchar(align)) {
    ch<-substr(align, i, i)
    if (ch=='C') {
      score<-score+match.score
      if (gap.opened) {gap.opened<-FALSE}
    } else if (ch=='M') {
      score<-score+mismatch.score
      if (gap.opened) {gap.opened<-FALSE}
    } else if ((ch=='I') | (ch=='D')) {
      score<-score+ifelse(gap.opened, gap.extend.score, gap.open.score)
      if (!gap.opened) {gap.opened<-TRUE}
    }
  }
  score
}

is.indel.chopped<-function(align) {
  chopped<-FALSE
  greg.results<-gregexpr('([ID]((-*)?))+', align)
  if (greg.results[[1]][[1]]>0) {
    first.indel.start<-greg.results[[1]][[1]]
    first.indel.len<-attr(greg.results[[1]], 'match.length')[[1]]
    
    last.indel.start<-greg.results[[1]][[length(greg.results[[1]])]]
    last.indel.len<-attr(greg.results[[1]], 'match.length')[[length(greg.results[[1]])]]
    
    first.indel<-substr(align, 0, first.indel.start+first.indel.len-1)
    last.indel<-substr(align, last.indel.start, nchar(align))
    
    chopped<-((score.align(first.indel)<0) | (score.align(last.indel)<0))
    #print(score.align(first.indel))
    #print(score.align(last.indel))
  }
  chopped
}
