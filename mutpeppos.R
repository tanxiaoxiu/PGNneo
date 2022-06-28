##Localization mutation sites of peptides
rm(list = ls(all = TRUE))
library(data.table)
library(stringr)

args <- commandArgs(T)
inFile <- args[1]
outFile <- args[2]
dt=fread(inFile)
dt=data.frame(dt)
rownames(dt) <- dt[,1]
dt[,1] <- NULL
dt['mut1']=NA
dt['mut11']=NA
dt['site1']=NA
dt['site11']=NA

dt['mut2']=NA
dt['mut22']=NA
dt['site2']=NA
dt['site22']=NA

dt['mut3']=NA
dt['mut33']=NA
dt['site3']=NA
dt['site33']=NA

for (i in 1:nrow(dt)) {

  if(nchar(dt[i,9])>0){
    key=33:(33+ceiling((nchar(dt[i,4])+2)/3)-1)
    mut=str_sub(dt[i,9],33,(33+ceiling((nchar(dt[i,4])+2)/3)-1))
    sp=data.frame(str_locate_all(dt[i,9],"\\*")[[1]])$start

    if (length(sp)>0) {

      xiao=which(sp<33)
      if (length(xiao)>0) {
        xiao=xiao[length(xiao)]
        start=sp[xiao]+1
        start2=sp[xiao]
      }else{
        start=1
        start2=0
      }

      da=which(sp>33)

      if(length(da)>0)
      {
        da=da[1]
        end=sp[da]-1
        end2=sp[da]
      }else{

        end=nchar(dt[i,9])
        end2=nchar(dt[i,9])

      }
    }else{
      start=1
      start2=0
      end=nchar(dt[i,9])
    }


    nseq=str_sub(dt[i,9],start,end)


    sp2=data.frame(str_locate_all(mut,"\\*")[[1]])$start

    yes=which(key[sp2]==33)
    if(length(yes)==0){


      if (length(sp2)>0) {

        dt[i,'mut1']=nseq
        dt[i,'mut11']=str_sub(dt[i,9],key[1],key[sp2-1][1])
        dt[i,'site1']=paste(start,end,sep='-')
        dt[i, 'site11']=paste((33-start2),( key[sp2-1][1]-start2),sep='-')

      }else{

        dt[i,'mut1']=nseq
        dt[i,'mut11']=mut
        dt[i,'site1']=paste(start,end,sep='-')
        dt[i, 'site11']=paste((33-start2),((33+ceiling((nchar(dt[i,4])+2)/3)-1)-start2),sep='-')
      }


    }



  }




  if(nchar(dt[i,8])>0){
    key=34:(34+ceiling(nchar(dt[i,4])/3)-1)
    mut=str_sub(dt[i,8],34,(34+ceiling(nchar(dt[i,4])/3)-1))
    sp=data.frame(str_locate_all(dt[i,8],"\\*")[[1]])$start

    if (length(sp)>0) {

      xiao=which(sp<34)
      if (length(xiao)>0) {
        xiao=xiao[length(xiao)]
        start=sp[xiao]+1
        start2=sp[xiao]
      }else{
        start=1
        start2=0
      }

      da=which(sp>34)

      if(length(da)>0)
      {
        da=da[1]
        end=sp[da]-1
        end2=sp[da]
      }else{

        end=nchar(dt[i,8])
        end2=nchar(dt[i,8])

      }
    }else{
      start=1
      start2=0
      end=nchar(dt[i,8])
    }




    nseq=str_sub(dt[i,8],start,end)


    sp2=data.frame(str_locate_all(mut,"\\*")[[1]])$start

    yes=which(key[sp2]==34)
    if(length(yes)==0){


      if (length(sp2)>0) {

        dt[i,'mut2']=nseq
        dt[i,'mut22']=str_sub(dt[i,8],key[1],key[sp2-1][1])
        dt[i,'site2']=paste(start,end,sep='-')
        dt[i, 'site22']=paste((34-start2),( key[sp2-1][1]-start2),sep='-')

      }else{

        dt[i,'mut2']=nseq
        dt[i,'mut22']=mut
        dt[i,'site2']=paste(start,end,sep='-')
        dt[i, 'site22']=paste((34-start2),((34+ceiling(nchar(dt[i,4])/3)-1)-start2),sep='-')
      }


    }



  }




  if(nchar(dt[i,7])>0){
    key=34:(34+ceiling((nchar(dt[i,4])+1)/3)-1)
    mut=str_sub(dt[i,7],34,(34+ceiling((nchar(dt[i,4])+1)/3)-1))
    sp=data.frame(str_locate_all(dt[i,7],"\\*")[[1]])$start

    if (length(sp)>0) {

      xiao=which(sp<34)
      if (length(xiao)>0) {
        xiao=xiao[length(xiao)]
        start=sp[xiao]+1
        start2=sp[xiao]
      }else{
        start=1
        start2=0
      }

      da=which(sp>34)

      if(length(da)>0)
      {
        da=da[1]
        end=sp[da]-1
        end2=sp[da]
      }else{

        end=nchar(dt[i,7])
        end2=nchar(dt[i,7])

      }
    }else{
      start=1
      start2=0
      end=nchar(dt[i,7])
    }




    nseq=str_sub(dt[i,7],start,end)


    sp2=data.frame(str_locate_all(mut,"\\*")[[1]])$start

    yes=which(key[sp2]==34)
    if(length(yes)==0){


      if (length(sp2)>0) {

        dt[i,'mut3']=nseq
        dt[i,'mut33']=str_sub(dt[i,7],key[1],key[sp2-1][1])
        dt[i,'site3']=paste(start,end,sep='-')
        dt[i, 'site33']=paste((34-start2),( key[sp2-1][1]-start2),sep='-')

      }else{

        dt[i,'mut3']=nseq
        dt[i,'mut33']=mut
        dt[i,'site3']=paste(start,end,sep='-')
        dt[i, 'site33']=paste((34-start2),((34+ceiling((nchar(dt[i,4])+1)/3)-1)-start2),sep='-')
      }


    }



  }



#  print(i)

}



colnames(dt)=c('chr','pso',"ref",'mut','gene','seq','prot1','prot2','prot3',
               'prot3-ant','prot3-mut','prot3-ant-pos','prot3-mut-pos',
               'prot2-ant','prot2-mut','prot2-ant-pos','prot2-mut-pos',
               'prot1-ant','prot1-mut','prot1-ant-pos','prot1-mut-pos')



write.csv(dt,file = outFile,quote = F)


