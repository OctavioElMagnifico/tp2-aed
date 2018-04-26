#Objetivo: Implementar el Descenso por el Gradiente

# Problemas a Resolver: 

# 1) Necesito algo afín al tipo ecuación. Parece que lo puedo solucionar con funciones que tomen dos vectores, uno de parámetros y uno de variables. 
# 2) De la derivación automática, mejor olvidarse por lo pronto. Para cada modelo voy a escribir la derivada. Sin embargo quizás con métodos numéricos hay un rebusque. 
# 3) Determinación del tamaño de paso. Para este, voy a implementar una búsqueda lineal por bipartición.  

library(tidyverse)

#cada modelo es un par de funciones que toma las mismas variables: nuestra función modeladora blahF con parámetros desconocidos y sus variables, y la derivada blahD de su error.



#Estasprimeras dos funciones tienen por finalidad representar a la familia de los modelos lineales y su derivada. 


mLineal <- function(x,parámetros)  {
    if (length(parámetros) != 2) {return("Cantidad equivocada de coficientes, este modelo toma 2.")}
    m  <- parámetros[1]
    b  <- parámetros[2]
    return(m*x+b)    
    }


gradienteML <-function(datos,parámetros){
    m   <- parámetros[1]
    b   <- parámetros[2]
    xs  <- datos[1]
    ys  <- datos[2]
    n   <- length(datos[[1]])
    fxs <- map(xs,~mLineal(.x,parámetros))
    variables <- c(xs,ys,fxs)
    dM  <- (2/n)  * sum( pmap_dbl(variables,~(..3-..2)*..1) )
    dB  <- (2/n)  * sum( pmap_dbl(variables,~(..3-..2)) )
    return(c(dM,dB))
    }

#Estas son funciones utilitarias para las partes del algoritmo

currificar <- function(modelo,parámetros,grad,h) {
    f <-function(x){
        u <- parámetros-h*grad
        return(modelo(x,u)) 
    }
    return(f)
}

normaEC <- function(x,y) {
    return((x-y)^2)
    }

errorF <- function(d,f) {
	x <- d[[1]]
	y <- d[[2]]
	n <- length(x)
	fxs <- map(x,f)
	v <- map2_dbl(fxs,y,normaEC)
	e <- sum(v)/n
	return(e)
}


#Estas son las partes del algoritmo

búsquedaPaso <- function(modelo,datos,parámetros,grad,cota) {
    exp <- -10
    e_0 <- 100
    e_1 <- 99
    paso <- 0
    while (e_1 < e_0) {
        paso=10^exp
        exp <- exp+1
        f <- currificar(modelo,parámetros,grad,10^exp)
        e_1 <- errorF(datos,f)
        }
    for (i in c(seq(0,paso,length=40))) {
        f <- currificar(modelo,parámetros,grad,i)
        e <- errorF(datos,f)
        if (e<e_1) {
            e_1 <- e
            paso <- i
            } 
        }
    return(paso)
}

descenderML <- function(datos,parámetros,cota) {
    gradiente <- c(0,0)
    e <- 100
    while (e>cota) {
        f <- currificar(mLineal,parámetros,gradiente,10)
        gradiente <- gradienteML(datos,parámetros)
        paso <- búsquedaPaso(mLineal,datos,parámetros,gradiente,cota)
        if (paso == 0) {return(parámetros)}
        parámetros <- parámetros - paso*gradiente
#        print(parámetros)
        e <- errorF(datos,f)
        }
    return(parámetros)
    }



autos <- read.table("/home/octavio/AnálisisExploratorio/Datos/autos.txt",TRUE)

tabla <- data.frame(x=autos$precio,y=autos$calidad)

descensoTabla <- descenderML(tabla,c(1,2),0.000001)
errorDescenso <- errorF(autos,~mLineal(.x,descensoTabla))
lmTabla <- unname(coefficients(lm(y~x+1,tabla)))
paramTabla <- c(lmTabla[2],lmTabla[1])
errorlm <- errorF(autos,~mLineal(.x,paramTabla))

