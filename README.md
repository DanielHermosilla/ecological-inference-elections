# Ecological Inference Elections

Esta versión es **provisoria** y provee funcionalidades limitadas en comparación a la versión final. En específico, todavía hace falta implementar el *"método bruto"* para la agregación de grupos. Las funciones `bootstrap()` y `get_group_agg()` deberían ser estables. `run_em()` también es estable, pero el log-likelihood presenta algunos errores. 

>[!INFO]
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

>[!ATTENTION]
> En caso de encontrar algún error, ayudaría **mucho** que me lo puedan hacer saber. Los *bugs* serán arreglados a la inmediatez y ayudarán con la versión final.
> Dejo mi contacto en caso de querer una atención más específica: **daniel.hermosilla.r@ug.uchile.cl**

## Uso


