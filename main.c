#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXIT   50  // Número máximo de iterações em métodos iterativos



typedef float real_t;

//Matriz tridiagonal
typedef struct {
    int n; //numero de pontos internos na malha
    double a, b; //intervalo
    double ya, yb;//condições de contorno
    double (* p)(double), (* q)(double), (* r)(double);
}Edo;

double getP(double x){
    return 0.0;
}

double getQ(double x){
    return 0.0;
}

double getR(double x){//Eq a.
    //6x - 0.5x²
    return 6*x - 0.5*x*x;
}

void gaussSeidel(Edo *edoeq, double *Y)
{
    int n = edoeq->n, k, i;
    double h, xi, bi, yi, d, di, ds;

    h = (edoeq->b - edoeq->a)/(n+1); //Largura do passo da malha
    for(k=0; k<50; ++k){
        for(i=0; i<n; ++i){
            xi = edoeq->a + (i+1)*h;//valor xi da malha
            bi = h*h * edoeq->r(xi);//termo independente
            di = 1 - h*edoeq->p(xi)/2.0;//diagonal inferior
            d = -2 + h*h * edoeq->q(xi);//diagonal principal
            ds = 1 + h*edoeq->p(xi)/2.0;//diagonal superior
            if(i==0)          bi -= ds*Y[i+1] + edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0);
            else if(i == n-1) bi -= di*Y[i-1] - edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0);
            else              bi -= ds*Y[i+1] + di*Y[i-1];

            Y[i] = bi/d; //calcula incognita

        }
    }

}




int main(){

    Edo *edoeq = (Edo *)malloc(sizeof(Edo));
    //y'' = 6x - 0.5x², x e (0,12)
    //y(0)=0 y(12) = 0
    edoeq->a = 0;
    edoeq->b = 12;
    edoeq->n = 5;
    edoeq->ya = 0;
    edoeq->yb = 0;
    edoeq->p = getP;
    edoeq->q = getQ;
    edoeq->r = getR;
    double *Y;
    Y = (double *)calloc((edoeq->n), sizeof(double));
    gaussSeidel(edoeq, Y);
    for(int i=0; i<edoeq->n; i++)
        printf("%f ", Y[i]);
    
    
}