# Ecological Inference Elections

Esta versión es **provisoria** y provee funcionalidades limitadas en comparación a la versión final. En específico, todavía hace falta implementar el *"método bruto"* para la agregación de grupos. Las funciones `bootstrap()` y `get_group_agg()` deberían ser estables. `run_em()` también es estable, pero el log-likelihood presenta algunos errores. 

>[!TIP]
> La documentación está disponible en un `fastei.pdf` en inglés. Favor verificar el documento para el uso de las funciones.

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

## Documentación

A continuación una documentación traducida y resumida de las funciones más importantes: 

1. `run_em()`: Corre el algoritmo para estimar las probabilidades agregadas:
- `object`: Objeto de tipo `eim` (creado con la función `eim()`) que tienen las matrices a calcular. Opcional en caso de otorgar las matrices `X` y `W` de forma manual o un archivo `json`.
- `X`: Matriz de candidatos, donde las filas corresponden a las urnas y las columnas a la indexación del candidato. Cada campo serían los votos que obtuvo cierto candidato para una urna en específico.
- `W`: Matriz de grupos, donde las filas corresponden a las urnas y las columnas a grupos demográficos (i.e partidos políticos, sexo, etc). Los campos son los votos que hizo cada grupo para cierta urna.
- `json_path`: Ruta con un archivo `.json` con la información de las votaciones. Este parámetro es opcional si se optó por usar alguna de las opciones anteriores (otorgar un objeto o las dos matrices).
- `method`: El método a utilizar para las estimaciones. Está disponible: `mult`, `mvn_cdf`, `mvn_pdf`, `hnr` y `exact`. Por defecto se utiliza `mult`, que es el más eficiente.
- `initial_prob`: Método para calcular la probabilidad inicial. Están disponibles: `group_proportional`, `proportional` y `uniform`. Por defecto se utiliza el primer método mencionado.
- `allow_mismatch`: Boolean que indica al algoritmo si hay un descuadre entre los votos de las mesas. Por defecto es `false`.
- `maxiter`: Máximas iteraciones permitidas, por defecto, 1000.
- `stop_threshold`: Valor indicador para terminar las iteraciones. Si la diferencia entre probabilidades es menor a `stop_threshold`, se deja de iterar. Por defecto es 0.01.
- `log_threshold`: Valor indicador para terminar las iteraciones, en base al log-likelihood. Si la diferencia entre log-likelihoods es menor a `log_threshold` se deja de iterar. En caso de utilizarse, el algoritmo deja de iterar en base a qué se cumple primero (`stop_threshold` o `log_threshold`). Por defecto viene desactivado, con un valor `-Inf`.
- `verbose`: Boolean que indica si hay que printear información que podría llegar a ser útil. Por defecto está desactivado.
- `step_size`: El tamaño del "paso" en caso de utilizar `method=hnr`, por defecto se utiliza 3000.
- `samples`: Cantidad de samples en caso de utilizar `method=hnr`, por defecto se utilizan 1000.
- `mc_method`: El método para simular la Normal Multivariada en `method=mvn_cdf`. Puede ser `genz` o `genz2`. Por defecto se utiliza `genz2`, que tiende a ser más rápida.
- `mc_error`: Threshold para detener la simulación de la Normal Multivariada, por defecto, 1e-6.
- `mc_samples`: Cantidad máxima de sampleos para la simulación, por defecto, 5000.

2. `bootstrap()`: Corre un algoritmo de bootstrapping para simular la desviación estándar. Hereda todos los argumentos de `run_em()` y adicionalmente considera:
- `nboot`: Cantidad de iteraciones a utilizar. Por defecto son 50.
- `seed`: La semilla para fijar las simulaciones. Por defecto, no existe impone ninguna semilla.

## Uso

A continuación se dará una breve demostración que lo pueden correr en su formato favorito (`.rmd`, `.ipynb`, `R`, `.qmd`, etc...). En caso de usar `.ipynb`, [instalar el Kernel de R](https://github.com/IRkernel/IRkernel).

Para este ejemplo utilizaremos los datos de la segunda vuelta de las elecciones presidenciales en Chile del $2021$. Afortunadamente, la librería ya viene con el dataset creado:

```r
library(fastei)

data("chile_election_2021")
election_data <- get("chile_election_2021") # Guardamos los datos en la variable
```

Si bien el dataset también está documentado en el `.pdf`, es posible ver lo que trae:

```r
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

Para este análisis, podría ser de interés estimar como vota el distrito electoral "El Golf":

```R
el_golf <- election_data[election_data$ELECTORAL.DISTRICT == "EL GOLF", ]

# Cantidad de urnas
nrow(el_golf) # 318
```

Bastaría dejar los resultados de los candidatos como una matrix de dimensión `(b x c)` (urnas como filas, candidatos como columnas) y los grupos como `(b x g)`. Luego, es directo pasar los datos a la función `run_em()`:

```R
candidatos <- as.matrix(el_golf[, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "BLANK.VOTES", "NULL.VOTES")])
grupos <- as.matrix(el_golf[, c("X18.19", "X20.29", "X30.39", "X40.49", "X50.59", "X60.69", "X70.79", "X80.")])

resultados <- run_em(X = candidatos, W = grupos)
# Error in run_em(X = candidatos, W = grupos, method = "mult") : 
#  run_em: Mismatch in the number of votes: 'X' has 67464 votes, but 'W' has 67450 votes. To allow a mismatch, set the argument to TRUE.
```

La función da un error. Esto se debe porque existe un error de descuadre entre los votos de los candidatos y grupos. En específico, las filas de ambas matrices no tienen las misma cantidad de votos.

Afortunadamente, es posible permitir correr el algoritmo de todas formas, pero debemos indicar explicitamente que lo estamos haciendo así. Lo anterior hace unos leves ajustes técnicos y nos restringe a utilizar únicamente los métodos `mult`, `mvn_cdf` y `mvn_pdf`.

```R
resultados <- run_em(X = candidatos, W = grupos, allow_mismatch = TRUE)
resultados

#   eim ecological inference model
#   Candidate matrix (X) [b x c]:
#   #   C1  C2 C3 C4 C5 C6 C7 BLANK.VOTES NULL.VOTES
#   32023 27 124 17 60  1  3  3           0          1
#   32024 28 129 12 74  0  3  0           0          3
#   32025 24 114 19 82  3  2  3           0          0
#   32026 24 126 16 65  1  3  3           1          1
#   32027 26 110 13 86  1  4  4           0          0
#   .
#   .
#   .
#   Group matrix (W) [b x g]:
#           X18.19 X20.29 X30.39 X40.49 X50.59 X60.69 X70.79 X80.
#   32023      8     17     72     50     34     29     19    7
#   32024      4     29    101     51     25     18     14    7
#   32025      3     23    100     48     25     31     11    6
#   32026      2     20     88     47     35     23     17    8
#   32027      3     15     93     60     29     24     14    6
#   .
#   .
#   .
#   Estimated probability [g x c]:
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#   [1,] 0.21 0.34 0.07 0.36 0.01 0.00 0.01 0.00 0.00
#   [2,] 0.23 0.34 0.05 0.34 0.00 0.01 0.01 0.01 0.01
#   [3,] 0.17 0.38 0.07 0.35 0.00 0.01 0.01 0.00 0.00
#   [4,] 0.08 0.49 0.08 0.33 0.00 0.01 0.00 0.00 0.00
#   [5,] 0.07 0.53 0.08 0.29 0.00 0.02 0.01 0.00 0.00
#   .
#   .
#   .
#   Method:  mult 
#   Total Iterations: 40 
#   Total Time (s): 0.012418 
#   Log-likelihood: 1.4788
```

Notemos que, al printear la matriz, sólo muestran algunos resultados, truncando las filas (y redondeando los valores al segundo decimal). Para acceder de forma exacta a cada atributo, se acceden como si fueran una lista:

```R
names(resultados)
# [1] "X"          "W"          "method"     "prob"       "logLik"     "iterations" "time"       "message"    "status"     "q"

resultados$prob
#           [,1]     [,2]      [,3]     [,4]        [,5]       [,6]        [,7]        [,8]        [,9]
#   [1,] 0.2082223 0.341486 0.0698584 0.358382 0.008261427 0.00322275 0.005523135 0.004337948 0.000706517
#   [2,] 0.2346563 0.335645 0.0503062 0.338688 0.004576561 0.00782078 0.010227206 0.008309655 0.009770271
#   [3,] 0.1708502 0.377271 0.0687785 0.352820 0.002913483 0.00913784 0.013066599 0.000849807 0.004312250
#   [4,] 0.0787467 0.493573 0.0782038 0.331677 0.001614309 0.00923532 0.004038647 0.000349432 0.002561771
#   [5,] 0.0717976 0.528725 0.0817829 0.291648 0.000448955 0.01649627 0.007733976 0.000625880 0.000741022
#   [6,] 0.0595029 0.562893 0.0613441 0.296641 0.002955823 0.01008786 0.004399579 0.000442334 0.001732946
#   [7,] 0.0415448 0.638584 0.0647837 0.236084 0.000919664 0.01371801 0.001014755 0.002345759 0.001006154
#   [8,] 0.0398162 0.642195 0.0660913 0.232982 0.005185203 0.00262540 0.000356071 0.002164199 0.008584368

resultados$iterations
# 40

resultados$message
# "Convergence achieved"
```

Supongamos que es deseable conocer una estimación a la desviación estándar. Es posible correr un algoritmo de [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)). Nos aseguramos de utilizar una semilla dado que la función hace muestreos aleatorios.

```R 
resultados <- bootstrap(resultados, nboot = 200, seed = 42)
resultados$sd

#           [,1]       [,2]       [,3]       [,4]        [,5]        [,6]        [,7]        [,8]        [,9]
#   [1,] 0.02407474 0.02724292 0.00753039 0.00976725 0.002978519 0.004003824 0.001695227 0.001714175 0.001397220
#   [2,] 0.00886291 0.01184461 0.00356447 0.00728163 0.000751909 0.001241132 0.001251894 0.001151976 0.001209765
#   [3,] 0.00700509 0.01159313 0.00345611 0.00699608 0.000702400 0.001266271 0.001215726 0.000492714 0.000807001
#   [4,] 0.00425882 0.00641300 0.00359046 0.00640130 0.000571504 0.001185845 0.000864958 0.000390415 0.000594928
#   [5,] 0.00401962 0.00660867 0.00301243 0.00502705 0.000327425 0.000978794 0.000697027 0.000527377 0.000505639
#   [6,] 0.00426301 0.00735525 0.00241111 0.00459194 0.000420586 0.000906268 0.000598853 0.000394278 0.000437642
#   [7,] 0.00574961 0.01221880 0.00324491 0.00657881 0.000442474 0.001353169 0.000685534 0.000475568 0.000489048
#   [8,] 0.00726865 0.01746242 0.00503719 0.00968632 0.000844869 0.001384006 0.000862021 0.000720071 0.000947733
```

Existen más funcionalidades anexas, como `write.csv(resultados, "resultados.json")` para guardar los resultados, o `summary()` para retornar una lista de los atributos más importantes.

