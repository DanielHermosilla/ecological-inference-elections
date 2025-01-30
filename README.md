# Estructura del repositorio 

El repositorio tiene la siguiente estructura:

```
.
├── DESCRIPTION
├── NAMESPACE
├── main.R
└── src
    ├── CMakeLists.txt
    ├── Makefile
    ├── analyzeResults.py
    ├── compile_commands.json
    ├── exact.c
    ├── exact.h
    ├── globals.h
    ├── hitAndRun.c
    ├── hitAndRun.h
    ├── main.c
    ├── main.h
    ├── matrixUtils.c
    ├── matrixUtils.h
    ├── multinomial.c
    ├── multinomial.h
    ├── multivariate-cdf.c
    ├── multivariate-cdf.h
    ├── multivariate-pdf.c
    ├── multivariate-pdf.h
    ├── tests
    │   ├── compile_commands.json
    │   └── test-matrixUtils.c
    └── utils
        ├── combinations.c
        ├── combinations.h
        ├── fileutils.c
        ├── fileutils.h
        ├── instanceGenerator.c
        ├── instanceGenerator.h
        ├── matrixUtils.c
        ├── matrixUtils.h
        ├── memoizationUtil.c
        ├── memoizationUtil.h
        ├── multivariateUtils.c
        ├── multivariateUtils.h
        └── uthash.h

4 directories, 37 files
```
Como forma generalizada, se tiene que: 

- `src`: Directorio con todo el código fuente a compilar. Ahí está todo el código `.c` y `.h`. Además, tiene un archivo `analyzeResults.py` para analizar los resultados.

- `instances`: Directorio con todas las instancias en formato `.json`. 

# Compilación

## Repositorio

Si es primera vez que utilizamos el repositorio, corremos: 

```
git clone https://github.com/DanielHermosilla/ecological-inference-elections.git 
```
Si queremos actualizarlo con los últimos cambios: 

```
git fetch --all 
```

Si además, se quieren borrar los cambios locales para tener un repositorio identico al de Github: 

```
git reset --hard origin/main 
``` 

## Preparación del directorio de binarios 

Aseguramos que `CMakeLists.txt` tenga las instrucciones locales de nuestro computador. De forma resumida, se deben tener las librerias gsl, openblas, cJSON, lapack y cblas (conocida por openblas). Se pueden descargar al correr: 

```
brew install gsl openblas cJSON lapack
```

Adicionalmente, recomiendo utilizar el compilador gcc: 

`brew install gcc`

Para verificar las rutas de las instalaciones se puede hacer mediante 

```
brew info openblas gsl cJSON lapack
```

Después de asegurar que tenemos todas las dependencias, creamos el directorio `build` dentro de `src` y entramos a él: 

```
cd src && mkdir build && cd build
```

Dentro de build se correran dos comandos; uno para definir el compilador y preparar los metadatos y otro para hacer la compilación en sí. Para el primero, revisar que se ponga la ruta del compilador gcc.

```
cmake .. -DCMAKE_C_COMPILER=/usr/local/Cellar/gcc/14.2.0_1/bin/gcc-14 -DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc/14.2.0_1/bin/g++-14
```

Con el siguiente comando se compila el proyecto:

```
cmake --build .
```

El ejecutable estará dentro del directorio build.

# Uso 

## Correr las instancias 

Al haber compilado, se tendrá algo del siguiente estilo: 

```
.
├── CMakeLists.txt
├── Makefile
├── analyzeResults.py
├── build
│   ├── exact.d
│   ├── exact.o
│   ├── hitAndRun.d
│   ├── hitAndRun.o
│   ├── libutil.dylib
│   ├── main.d
│   ├── main.o
│   ├── multinomial.d
│   ├── multinomial.o
│   ├── multivariate-cdf.d
│   ├── multivariate-cdf.o
│   ├── multivariate-pdf.d
│   ├── multivariate-pdf.o
│   ├── util_exec
│   └── utils
│       ├── combinations.d
│       ├── combinations.o
│       ├── fileutils.d
│       ├── fileutils.o
│       ├── instanceGenerator.d
│       ├── instanceGenerator.o
│       ├── matrixUtils.d
│       ├── matrixUtils.o
│       ├── memoizationUtil.d
│       ├── memoizationUtil.o
│       ├── multivariateUtils.d
│       └── multivariateUtils.o
├── compile_commands.json
├── exact.c
├── exact.h
├── globals.h
├── hitAndRun.c
├── hitAndRun.h
├── main.c
├── main.h
├── matrixUtils.c
├── matrixUtils.h
├── multinomial.c
├── multinomial.h
├── multivariate-cdf.c
├── multivariate-cdf.h
├── multivariate-pdf.c
├── multivariate-pdf.h
├── tests
│   ├── compile_commands.json
│   └── test-matrixUtils.c
└── utils
    ├── combinations.c
    ├── combinations.h
    ├── fileutils.c
    ├── fileutils.h
    ├── instanceGenerator.c
    ├── instanceGenerator.h
    ├── matrixUtils.c
    ├── matrixUtils.h
    ├── memoizationUtil.c
    ├── memoizationUtil.h
    ├── multivariateUtils.c
    ├── multivariateUtils.h
    └── uthash.h

5 directories, 60 files
```

Para correr el ejecutable, hay que estar dentro del directorio `build`.  La estructura del llamado es la siguiente: 

```
./util_exec "[directorio de instancias]" "[directorio de resultados]" "[Método]"
```

A modo de ejemplo, si quiero correr el ejecutable con el método multinomial, guardar los resultados en una carpeta llamada results y **utilizar las instancias que ya vienen dentro del repositorio**, es posible correr: 

```
./util_exec "../../instances" "results" "Multinomial"
```
Lo que guardará los resultados en `/build/results/Multinomial`. En caso de haber un error, asegurar crear el directorio de resultados con `mkdir`. 

## Análisis de resultados 

Para poder analizar los resultados de forma agrupada por grupo y candidato, se provee el script bajo el nombre `src/analyzeResults.py`. Su uso es el siguiente: 

```
python analyzeResults.py "[directorio de los resultados de instancias]" "[directorio con los resultados analizados]" "[Método]"
```

El directorio con los resultados de las instancias es el mismo al que pasamos como segundo argumento al llamar el ejecutable. 

A modo de ejemplo, ocupando el mismo comando del inciso anterior: 

```
python analyzeResults.py "build/results" "finalResults" "Multinomial"
```
Creará una carpeta en `project/src/finalResults` con los resultados promediados. 
