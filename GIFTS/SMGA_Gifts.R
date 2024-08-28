###Script to pars the DRAM annotation and calculate GIFTs and MCI

#Author: Arturo Vera

library(tidyverse)
library(distillR)
library(Rtsne)
library(vegan)
library(patchwork)
library(viridis)
library(Polychrome)
library(broom)

#Load the Annotation table



SMGAAnnot <- readRDS("SMGA_Total_Annot_Metadata.RDS") 


SMGAtoGIFTs <- SMGAAnnot %>%
  select(Genome,Gene_id,kegg_id)


GIFTs <- distill(SMGAtoGIFTs,GIFT_db,genomecol=1,annotcol=3)


#Aggregate bundle-level GIFTs into the compound level
GIFTs_elements <- to.elements(GIFTs,GIFT_db)

#Aggregate element-level GIFTs into the function level
GIFTs_functions <- to.functions(GIFTs_elements,GIFT_db)

#Aggregate function-level GIFTs into overall Biosynthesis, Degradation and Structural GIFTs
GIFTs_domains <- to.domains(GIFTs_functions,GIFT_db)

avg_mci <- rowMeans(GIFTs_functions) %>%
  enframe(name = "Genome",value = "MCI")

#Get Tax info and water source


Metadata <- SMGAAnnot %>%
  select(Genome,Order,Genus,WaterSource,Color) %>%
  distinct(Genome,.keep_all = T)

#Join the tables

avg_mci_Meta <- avg_mci %>%
  left_join(Metadata,by="Genome")


##tsne
set.seed(123)

tsne_func <- Rtsne(X = GIFTs_elements, 
                   dims = 2, 
                   check_duplicates = FALSE)


tsne_df <- as_tibble(tsne_func$Y) %>%
  rename(tSNE1 = V1, tSNE2 = V2) %>%
  bind_cols(avg_mci_Meta) %>%
  select(Genome,everything())
  

TSNEplotOrder <- tsne_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2, 
             color = Order,
             shape=WaterSource)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values=(tsne_df %>%
                                        select(Order,Color) %>%
                                        distinct() %>%
                                        deframe())) +
  theme_minimal()

##Add The genus labels

SUB <- tsne_df %>%
  group_by(Genus) %>%
  summarize(Genome=first(Genome)) %>%
  select(Genome) %>%
  left_join(tsne_df,by="Genome")

TSNEplotOrderGenus <- tsne_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2, 
             color = Order,
             shape=WaterSource)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values=(tsne_df %>%
                               select(Order,Color) %>%
                               distinct() %>%
                               deframe())) +
  geom_text(data=SUB,
            aes(x=tSNE1,
                y=tSNE2,
                label=Genus),
            color="black",
            show.legend = F,
            hjust = 1,
            vjust=2) +
  geom_segment(data = SUB, 
               aes(x = tSNE1,
                   y = tSNE2,
                   xend = tSNE1,
                   yend = tSNE2 - 1), 
               color = "black", 
               size = 0.5)+
  theme_minimal()


TSNEplotMCI <- tsne_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2, 
             color = MCI,
             shape=WaterSource)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_gradientn(colours =  viridis(18)) +
  theme_minimal()


TSNECombined <- TSNEplotOrderGenus / TSNEplotMCI


ggsave(TSNECombined,
       file="TSNE.SMGA.pdf",
       width = 15,
       height = 15)
####Get the heatmap

Elements <- GIFT_db %>%
  select(Code_element,Domain, Function)

GIFTs_elements %>%
  pheatmap::pheatmap(cluster_cols = F)

#Add taxonomy to this table

GIFTs_W_Tax <- as_tibble(GIFTs_elements,rownames = NA) %>% 
  rownames_to_column("Genome")  %>%
  left_join(Metadata,by="Genome")

GIFTSannotCol <- Elements %>%
  select(Code_element,Function) %>%
  distinct() %>%
  column_to_rownames("Code_element")


GIFTSTaxAnnot <- GIFTs_W_Tax %>%
  select(Genome,Order) %>%
  column_to_rownames("Genome")

GIFTSWaterAnnot <- GIFTs_W_Tax %>%
  select(Genome,WaterSource) %>%
  column_to_rownames("Genome")
  

ColorAnnot <- list(Order=(Metadata %>%
                            select(Order,Color) %>%
                            distinct() %>%
                            deframe()),
                   WaterSource=c("blue","red"),
                   Function=(Elements %>%
                               select(Function) %>%
                               distinct() %>%
                               mutate(Color=palette36.colors(n=21)) %>%
                               deframe()))


Color <- rev(viridis(10))
GIFTsElementsPH <- GIFTs_elements %>%
  pheatmap::pheatmap(cluster_cols = F,
                     color = Color,
                     annotation_row = GIFTSTaxAnnot,
                     annotation_col = GIFTSannotCol,
                     annotation_colors = ColorAnnot,
                     show_colnames = F,
                     show_rownames = F)
ggsave(GIFTsElementsPH,
       file="GIFTsElementsPH.pdf",
       width = 12,
       height = 12)

#Calculate PERMANOVA

Perma <- tsne_df %>%
    select(Genome,MCI,WaterSource) %>%
    column_to_rownames("Genome")


PermaMat <- Perma %>%
  select(MCI)

permanova_result <- adonis2(PermaMat ~ Perma$WaterSource, method = "euclidean")

broom::tidy(permanova_result)

summary_df <- as.data.frame(summary(permanova_result))
