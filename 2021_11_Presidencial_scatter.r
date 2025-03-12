rm(list=ls())
cat("\014")

umbral = 1E-5
presicion = 1E-9

setwd("C:/Users/pablo/OneDrive/Desktop/TESIS/Tesis/model_v2")
PATH_READ_FILES = paste0(getwd(), "/2021_11_Presidencial/output/")
PATH_WRITE_FILES = paste0(getwd(), "/images/elections/")

library(ggplot2)
library(scales)
library(wesanderson)
require("ggrepel")
set.seed(42)
df_pais <- read.csv(paste0(PATH_READ_FILES, '2021_11_Presidencial_PAIS.csv'), sep = ',', header = TRUE, encoding = "UTF-8")
df_pais$CE_mesa = paste0(df_pais$CIRCUNSCRIPCION.ELECTORAL, "-", df_pais$MESA)
df_pais$P.VALOR = pmax(df_pais$P.VALOR, presicion/10)
# df_pais$DESCUADRADA = as.factor(df_pais$DESCUADRADA)
# df_pais$DESCUADRADA <- relevel(df_pais$DESCUADRADA, "1")
# levels(df_pais$DESCUADRADA) = c("Sí", "No")
# summary(df_pais)

if(F) {
  # grafico para presetnación del SERVEL
  ggplot(df_pais, aes(x=NUM.VOTOS, y=P.VALOR, color = NUM.MESAS)) + 
    geom_point(size = 1) + theme_dark() +
    xlab("Número de votos") + ylab("p-valor") +
    scale_y_log10() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(2*presicion/100, 1))  +
    scale_color_gradientn(colours = wes_palette("Zissou1", 10, type = "continuous"), name = "Número de mesas en local") + 
    theme(legend.position="top") + 
    guides(colour=guide_colourbar(barwidth=15,label.position="top")) + 
    geom_text_repel(data = subset(df_pais, P.VALOR < umbral), aes(label = CE_mesa),
                    size = 2, max.overlaps = 1000, fontface = "bold") 
  
  ggsave(paste0(PATH_WRITE_FILES, "scatter_pval_log_col_w.pdf"), width = 1800, height = 2450, units = "px")
}
#number_text = c("1 in ten", "1 in one hundred", "1 en a thousand", "1 in ten thousand", 
#        "1 in a hundred thousand", "1 in a million", "1 in ten millions", "1 a hundred million", 
#         "1 in a billion", "1 en ten billions", "< 1 in ten billions")
#
number_text = c("1 in ten", "1 in one hundred", "1 in one thousand", "1 in ten thousand", 
                "1 in one hundred \nthousand", "1 in one million", "1 in ten million", "1 in one hundred \nmillion",
                "1 in one billion", "< 1 in one billion")

#number_text = c("1 in ten", "1 in one hundred", "1 in one thousand", "1 in ten thousand", 
#                "1 in one hundred \nthousand", "1 in one million", "1 in ten million", "< 1 in ten million")
df_texto = data.frame(x = rep(485, -log10(presicion)+1), y = 10^seq(-1,log10(presicion)-1), 
                      label = number_text)


#colores = colorRampPalette(c("#BEEFAF", "#FFFB9F", "#ff9797"))(nrow(df_texto))
colores = colorRampPalette(c("#BEEFAF", "#FFFB9F", "#FFBB77", "#ff9797", "#EC8EA7"))(nrow(df_texto))
  

# grafico para presetnación del SERVEL
ggplot(df_pais, aes(x=NUM.VOTOS, y=P.VALOR)) + 
    geom_point(size = 1.5, color = 'blue3', alpha = 0.28, shape = 19, stroke=1) + theme_bw() +
    xlab("Number of votes") + ylab("p-value") +
    scale_y_log10(breaks = 10^seq(0,log10(presicion)-1),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(8*presicion/100, 1))  +
    #scale_color_manual(values = wes_palette("Darjeeling1"), name = "Mesa descuadrada") + 
    theme(legend.position="top", legend.background = element_rect(fill = "white", color = "black")) + 
    coord_cartesian(xlim = c(0, 450), # This focuses the x-axis on the range of interest
                 clip = 'off') +   # This keeps the labels from disappearing
    theme(plot.margin = unit(c(1,8,1,0), "lines")) + # This widens the right margin
    geom_label(data=df_texto, aes( x=x, y=y, label=label),
             color = "black", size=3, fill = colores, 
             hjust = "left") +
    geom_hline(yintercept=df_texto$y, linetype="dashed", color = "#404040")
    # geom_text_repel(data = subset(df_pais, P.VALOR < umbral/1000000), aes(label = CE_mesa),
    #                 size = 2, max.overlaps = 1000, fontface = "bold") 

#ggsave(paste0(PATH_WRITE_FILES, "scatter_pval_log_col_w_total.pdf"), width = 2000, height = 2450, units = "px")
#ggsave(paste0(PATH_WRITE_FILES, "scatter_pval_log_col_w_total.png"), width = 2000, height = 2450, units = "px")




