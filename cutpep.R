##get 8-11mer peptides
args<-commandArgs(T)
setwd(args[1])
library(stringr)
dat<-read.csv('maxqpep.txt',sep='\t')
for(len in 8:11)
{
  a<-data.frame(ff="")
  f=1
  name_end<-1
  for(k in 1:nrow(dat))
  {
    Position<-dat[k,'Mut.Position']
    len_seq<-str_length(dat[k,'Sequence'])
    start<-Position-len+1
    end<-Position
    for(i in 1:len)
    {
      if(len_seq>=len)
      {
        if(start<=0)
        {
          end=end+abs(start)+1
          start=1
          #print(str_sub(dat[k,'Sequence'],start,end))
          a[f,1]<-str_c('>',dat[k,2],'|',dat[k,3],'|',dat[k,6],'|',name_end)
          a[f+1,1]<-str_sub(dat[k,'Sequence'],start,end)
          f=f+2
          name_end=name_end+1
          #print(2)
          start=start+1
          end=end+1
        }
        else
        {
          #print(str_sub(dat[k,'Sequence'],start,end))
          a[f,1]<-str_c('>',dat[k,2],'|',dat[k,3],'|',dat[k,6],'|',name_end)
          a[f+1,1]<-str_sub(dat[k,'Sequence'],start,end)
          f=f+2
          name_end=name_end+1
          start=start+1
          end=end+1
        }
        if(end>len_seq)
        {
          break
        }
        if(start>Position)
        {
          break
        }
      }
    }
  }
  print(name_end)
  write.table(a,str_c('maxqpep_',len,'.fasta'),quote = F,row.names = F,col.names = F)
}
