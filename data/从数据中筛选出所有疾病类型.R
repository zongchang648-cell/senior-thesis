library(readxl)
library(writexl)
the_data<-read_xlsx('C:\\毕业论文\\数据\\PTEN\\uniprot_processed.xlsx',sheet=2)
total<-paste(the_data$disease,collapse = "/")
the_data2<-strsplit(total,"/")[[1]]
the_data2<-unique(the_data2)
the_data2<-as.data.frame(the_data2)
write_xlsx(the_data2,'C:\\毕业论文\\数据\\PTEN\\disease_type.xlsx')