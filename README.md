# Ecological Inference Elections

Esta versión es **provisoria** y provee funcionalidades limitadas en comparación a la versión final. En específico, todavía hace falta implementar el *"método bruto"* para la agregación de grupos. Las funciones `bootstrap()` y `get_group_agg()` deberían ser estables. `run_em()` también es estable, pero el log-likelihood presenta algunos errores. 

>[!TIP]
> La documentación está disponible en un `.pdf` en Inglés.

## Instalación

Gran parte del lenguaje está hecho por detrás en C. En realidad, R en sí está hecho en C. Por lo tanto, no se está *"reinventando la rueda"* al instalar este paquete.

R posee una forma nativa para instalar la librería de forma directa en Github. En teoría, bastaría con poner en la consola de R:

```{r}
devtools::install_github("DanielHermosilla/ecological-inference-elections")
```

Para poder trabajar con la librería, hay que importarla con `library(fastei)`. 

## Errores conocidos *"(Known bugs)"*

Dado la forma en que C administra la memoría, algunas veces existen crasheos al abortar en la mitad de generación de sampleos. **Esto tiene solución**, pero falta implementarla. Esto se va a arreglar, pero no está primero en la lista de prioridades. Por lo mismo, se desactivó, provisoriamente, el paralelismo en la generación de samples.

>[!WARNING]
> En caso de encontrar algún error, ayudaría **mucho** que me lo puedan hacer saber. Los *bugs* serán arreglados a la inmediatez y ayudarán con la versión final.
> Dejo mi contacto en caso de querer una atención más específica: **daniel.hermosilla.r@ug.uchile.cl**

## Uso

A continuación se dará una breve demostración que lo pueden correr en su formato favorito (`.rmd`, `.ipynb`, `R`, `.qmd`, etc...). En caso de usar `.ipynb`, [instalar el Kernel de R](https://github.com/IRkernel/IRkernel).

Para este ejemplo utilizaremos los datos de la segunda vuelta de las elecciones presidenciales en Chile del $2021$. Afortunadamente, la librería ya viene con el dataset creado:

```{r}
library(fastei)

data("chile_election_2021")
election_data <- get("chile_election_2021") # Guardamos los datos en la variable
```

Si bien el dataset también está documentado en el `.pdf`, es posible ver lo que trae:

```{r}
head(election_data)

#       ELECTORAL.DISTRICT BALLOT.BOX C1 C2 C3 C4 C5 C6 C7 BLANK.VOTES NULL.VOTES X18.19 X20.29 X30.39 X40.49 X50.59 X60.69 X70.79 X80. MISMATCH
#   1  ANTOFAGASTA NORTE        126 12 22  3  4  1  5 71           0          0      1      7    103      3      0      3      1    0    FALSE
#   2  ANTOFAGASTA NORTE        127 11 25  4  4  0  6 70           1          0      2     10     98      6      4      1      0    0     TRUE
#   3  ANTOFAGASTA NORTE        128 16 23  3  3  1  5 70           1          2      3      5    111      2      2      0      1    0     TRUE
#   4  ANTOFAGASTA NORTE        129 10 22  4  4  0  2 75           0          2      2      2    103      7      3      2      0    0     TRUE
#   5  ANTOFAGASTA NORTE        130 16 23  4  3  1  3 76           0          1      4      9    112      1      2      0      0    0     TRUE
#   6  ANTOFAGASTA NORTE        131 22 25  3  4  3  2 62           0          3      3      8    103      4      3      1      0    2     TRUE

summary(election_data)

#   ELECTORAL.DISTRICT  BALLOT.BOX              C1              C2              C3             C4              C5              C6             C7         BLANK.VOTES       NULL.VOTES        X18.19      
#   Length:46606       Length:46606       Min.   :  0.0   Min.   :  0.0   Min.   : 0.0   Min.   :  0.0   Min.   : 0.00   Min.   : 0.0   Min.   :  0.0   Min.   : 0.000   Min.   : 0.00   Min.   :  0.00  
#   Class :character   Class :character   1st Qu.: 25.0   1st Qu.: 29.0   1st Qu.:12.0   1st Qu.: 12.0   1st Qu.: 1.00   1st Qu.: 7.0   1st Qu.: 10.0   1st Qu.: 0.000   1st Qu.: 0.00   1st Qu.:  2.00  
#   Mode  :character   Mode  :character   Median : 36.0   Median : 39.0   Median :17.0   Median : 17.0   Median : 2.00   Median :11.0   Median : 17.0   Median : 0.000   Median : 1.00   Median :  3.00  
                                    #   Mean   : 38.5   Mean   : 41.9   Mean   :17.4   Mean   : 19.2   Mean   : 2.19   Mean   :11.4   Mean   : 19.3   Mean   : 0.652   Mean   : 1.19   Mean   :  4.24  
                                    #   3rd Qu.: 49.0   3rd Qu.: 50.0   3rd Qu.:22.0   3rd Qu.: 23.0   3rd Qu.: 3.00   3rd Qu.:15.0   3rd Qu.: 25.0   3rd Qu.: 1.000   3rd Qu.: 2.00   3rd Qu.:  5.00  
                                    #   Max.   :171.0   Max.   :200.0   Max.   :95.0   Max.   :118.0   Max.   :20.00   Max.   :49.0   Max.   :100.0   Max.   :12.000   Max.   :95.00   Max.   :246.00  
#
#   X20.29          X30.39          X40.49          X50.59          X60.69          X70.79           X80.     MISMATCH      
#   Min.   :  0.0   Min.   :  0.0   Min.   :  0.0   Min.   :  0.0   Min.   :  0.0   Min.   :  0.0   Min.   : 0   Mode :logical  
#   1st Qu.: 11.0   1st Qu.:  6.0   1st Qu.:  6.0   1st Qu.: 19.0   1st Qu.: 14.0   1st Qu.:  6.0   1st Qu.: 1   FALSE:9088     
#   Median : 16.0   Median : 20.0   Median : 23.0   Median : 28.0   Median : 24.0   Median : 12.0   Median : 3   TRUE :37518    
#   Mean   : 25.1   Mean   : 27.9   Mean   : 26.2   Mean   : 28.2   Mean   : 23.5   Mean   : 12.7   Mean   : 4                  
#   3rd Qu.: 26.0   3rd Qu.: 43.0   3rd Qu.: 41.0   3rd Qu.: 37.0   3rd Qu.: 32.0   3rd Qu.: 18.0   3rd Qu.: 6                  
#   Max.   :237.0   Max.   :166.0   Max.   :167.0   Max.   :182.0   Max.   :153.0   Max.   :125.0   Max.   :51                  
```


