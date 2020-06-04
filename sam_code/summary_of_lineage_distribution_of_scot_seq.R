##2020-05-22 summary of lineage distribution of scottish sequences
# lu lu
###note!!
#the number of global lineage from COG metadata table is different from Pangolin website, here use the data from COG table only.
#######################
#######################
library('dplyr')
library('scales')

path4<-'/Users/lulu/Documents/work/COVID19/phylodynamics of COVID19/reference data for lineage classcification/UK/'
name4<-'cog_global_2020-05-15_metadata'
FileName4 <- paste(path4, name4, ".csv", sep="")
ori_tbl1 <- read.csv(FileName4) 
ori_tbl1=as.data.frame(ori_tbl1)

#unique(ori_tbl1$country)

ori_tbl1$subcountry<-apply(as.matrix(ori_tbl1$sequence_name), 1, getEl, ind=1, sep="/")

glo_l<-ori_tbl1 %>% group_by(subcountry,country,lineage) %>% 
  summarize(freq = length(sequence_name)) %>% 
  ungroup()  %>%
  group_by(lineage) %>% summarise(all_freq=sum(freq))

UK_l<-ori_tbl1 %>% group_by(subcountry,country,lineage) %>% filter(country=='UK') %>% 
  summarize(freq = length(sequence_name)) %>% 
  ungroup()  %>%
  group_by(lineage) %>% summarise(UK_freq=sum(freq))

Scot_l<-ori_tbl1 %>% group_by(subcountry,country,lineage) %>% filter(subcountry=='Scotland') %>% 
  summarize(freq = length(sequence_name)) %>% 
  ungroup()  %>%
  group_by(lineage) %>% summarise(Scot_freq=sum(freq))


lineage_tbl<-glo_l %>% left_join(., UK_l) %>%
  left_join(.,Scot_l) 

write.csv(lineage_tbl, paste(path4,'lineage_distribution_global_UK_scot.csv'))

# scotland lineage
lineage_tbl$Scot_freq %>% is.na %>% `!` %>% sum()
scot <- lineage_tbl %>% dplyr::filter(!is.na(Scot_freq))%>% filter(!lineage %in% c('A','B'))
scot_lineage_count<-nrow(scot)
scot_lineages <- scot$lineage

##unique scot lineage
scot_lineage_uniq <- lineage_tbl %>% filter(UK_freq==Scot_freq) 


uk <- lineage_tbl %>% dplyr::filter(!is.na(UK_freq))%>% filter(!lineage %in% c('A','B'))
uk_lineage_count<-nrow(uk)
uk_lineages <- uk$lineage

####global lineage
glo <- lineage_tbl %>% dplyr::filter(!is.na(all_freq))%>% filter(!lineage %in% c('A','B'))
all_lineage_count<-nrow(glo)
all_lineages <- glo$lineage

######fraction
fraction_scot_uk<-percent(scot_lineage_count/uk_lineage_count)
fraction_scot_all<-percent(scot_lineage_count/all_lineage_count)


######print the summary
sprintf("Number of lineages in scotland = %s (They are %s). Fraction those are of all sublineages is %s; fraction those are of UK sublineages is %s. There is %s unique lineage to the UK, which is %s.", 
        scot_lineage_count, 
        scot_lineages%>%str_c(collapse = ','),
        fraction_scot_all,
        fraction_scot_uk,
        nrow(scot_lineage_uniq),
        scot_lineage_uniq$lineage) %>% cat()

