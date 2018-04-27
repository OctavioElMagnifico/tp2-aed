#Objetivo: Implementar el Descenso por el Gradiente

# Problemas a Resolver: 

# 1) Necesito algo afin al tipo ecuacion. Parece que lo puedo solucionar con funciones que tomen dos vectores, uno de parametros y uno de variables. 
# 2) De la derivacion automatica, mejor olvidarse por lo pronto. Para cada modelo voy a escribir la derivada. Sin embargo quizas con metodos numericos hay un rebusque. 
# 3) Determinacion del tamaño de paso. Para este, voy a implementar una busqueda lineal por biparticion.  

library(tidyverse)

#cada modelo es un par de funciones que toma las mismas variables: nuestra funcion modeladora blahF con parametros desconocidos y sus variables, y la derivada blahD de su error.



#Estasprimeras dos funciones tienen por finalidad representar a la familia de los modelos lineales y su derivada. 


mLineal <- function(x,parametros)  {
    if (length(parametros) != 2) {return("Cantidad equivocada de coficientes, este modelo toma 2.")}
    m  <- parametros[1]
    b  <- parametros[2]
    return(m*x+b)    
    }


gradienteML <-function(datos,parametros){
    m   <- parametros[1]
    b   <- parametros[2]
    xs  <- datos[1]
    ys  <- datos[2]
    n   <- length(datos[[1]])
    fxs <- map(xs,~mLineal(.x,parametros))
    variables <- c(xs,ys,fxs)
    dM  <- (2/n)  * sum( pmap_dbl(variables,~(..3-..2)*..1) )
    dB  <- (2/n)  * sum( pmap_dbl(variables,~(..3-..2)) )
    return(c(dM,dB))
    }

#Estas son funciones utilitarias para las partes del algoritmo

currificar <- function(modelo,parametros,grad,h) {
    f <-function(x){
        u <- parametros-h*grad
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

#aparetemente, el parametro con mayor influencia sobre la calidad del Descenso es la cantidad de partes en que el for parte el intervalo para buscar el paso, una vez hallada la cota superior.
#Esta cantidad departes claramente limita la precision del Descenso, sin importar el error elegido (en realidad se llega al return correspondiente a la repeticion del tamaño de pas0 0)

busquedaPaso <- function(modelo,datos,parametros,grad,cota,e_0) {
    exp <- -12
    e_1 <-  0
    paso <- 0
    while (e_1 < e_0) {
        paso=10^exp
        exp <- exp+1
        f <- currificar(modelo,parametros,grad,10^exp)
        e_1 <- errorF(datos,f)
        }
    for (i in c(seq(0,paso,length=500))) {
        f <- currificar(modelo,parametros,grad,i)
        e <- errorF(datos,f)
        if (e<e_1) {
            e_1 <- e
            paso <- i
            } 
        }
    return(c(paso,e_1))
}

descenderML <- function(datos,parametros,cota) {
    gradiente <- c(0,0)
    e <- 100
    while (e>cota) {
        gradiente <- gradienteML(datos,parametros)
        pasoyErr <- busquedaPaso(mLineal,datos,parametros,gradiente,cota,e)
        paso <- pasoyErr[1]
        ePrevio <- e
        e <- pasoyErr[2]
#        print("Error")
#        print(c(e,ePrevio))
        parametros <- parametros - paso*gradiente
        if ( abs(ePrevio - e) < cota ) {return(parametros)}
#        print("parametros")
#        print(parametros)
        }
    return(parametros)
    }



autos <- read.table("./autos.txt",TRUE)

tabla <- data.frame(x=autos$precio,y=autos$calidad)

mediaPrecio <- mean(autos$precio)
mediaCalidad <- mean(autos$calidad)

pendMedias <- (mediaCalidad / mediaPrecio)
ordMin <- min(autos$calidad)
vec0 <- c(ordMin,pendMedias)

descensoTabla <- descenderML(tabla,c(1,2),10^-7)
errorDescenso <- errorF(tabla,~mLineal(.x,descensoTabla))
lmTabla <- unname(coefficients(lm(y~x+1,tabla)))
paramTabla <- c(lmTabla[2],lmTabla[1])
errorlm <- errorF(autos,~mLineal(.x,paramTabla))
print("El error para autos del modelo del Descenso es")
print(errorDescenso)
print("Y el de la funcion incorporada")
print(errorlm)

archivoGalton <- read.csv("./GaltonMod.csv",TRUE)
galtonAtip <- data.frame(padre=archivoGalton$parent,hijo=archivoGalton$child)
galton <- galtonAtip[-4,]

descensoGalton <- descenderML(galton,c(1,28),10^(-7))
errGaltDesc <- errorF(galton,~mLineal(.x,descensoGalton))
lmGalton <- unname(coefficients(lm(hijo~padre+1,galton)))
parLmGalton <- c(lmGalton[2],lmGalton[1])
errorGaltLm <- errorF(galton,~mLineal(.x,parLmGalton))
print(c(errGaltDesc,errorGaltLm))
print("Arriba error descenso y error lm, abajo parámetros descenso y parámetros lm")
print(c(descensoGalton,parLmGalton))
