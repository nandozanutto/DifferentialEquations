#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXIT   50  // Número máximo de iterações em métodos iterativos

//y'' = 6x - 0.5x², x e (0,12)
//y(0)=0 y(12) = 0

typedef float real_t;

//Matriz tridiagonal
typedef struct {
    int n; //numero de pontos internos na malha
    double a, b; //intervalo
    double ya, yb;//condições de contorno
    double (* p)(double), (* q)(double), (* r)(double);
}Edo;


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
    int n=5;
    int a=0, b=12;
    real_t malha[n+2];
    real_t h;
    h = (b-a)/(n+1);
    malha[0] = a;
    malha[n+1] = b;
    for(int i=1; i<n+1; i++)
        malha[i] = a + i*h;
    
    // for(int i=0; i<n+2; i++)
    //     printf("%f\n", malha[i]);



}