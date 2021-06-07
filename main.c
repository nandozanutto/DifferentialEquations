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

typedef struct {
    int n, m; //numero de pontos internos na malha
    double Lx, Ly; //intervalo
    double (* u1)(double), (* u2)(double), (* u3)(double), (* u4)(double); //condições de contorno
    double (* func)(double, double);
}Edo2;

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

double func(double x, double y){//Eq a.
    //6x - 0.5x²
    return 0;
}

double u1(double x, double y){
    return 0;
}

double u2(double x, double y){
    return 0;
}

double u3(double x, double y){
    return 0;
}

double u4(double x, double y){
    return 0;
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
            else if(i == n-1) bi -= di*Y[i-1] + edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0);//Sinal - ou + ?
            else              bi -= ds*Y[i+1] + di*Y[i-1];

            Y[i] = bi/d; //calcula incognita

        }
    }

}

void gaussSeidel2(Edo2 *edoeq, double *Y)
{
    int n = edoeq->n, m = edoeq->m, i, k, j;
    double hx, hy, xi, bi, yi, d, di, ds, di2, ds2;

    hx = edoeq->Lx/(n+1);
    hy = edoeq->Ly/(m+1);//largura do passo
    
    for(k=0; k<50; ++k){
        
        for(j=0; j<m; j++){
            
            for(i=0; i<n; ++i){
                xi = (i+1)*hx;//valor xi da malha
                yi = (j+1)*hy;//valor yi da malha

                bi = hx*hx*hy*hy*edoeq->func(xi, yi);//termo indep.
                di = 1 - h*edoeq->p(xi)/2.0;//diagonal inferior
                d = -2 + h*h * edoeq->q(xi);//diagonal principal
                ds = 1 + h*edoeq->p(xi)/2.0;//diagonal superior
                if(i==0)          bi -= ds*Y[i+1] + edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0);
                else if(i == n-1) bi -= di*Y[i-1] + edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0);//Sinal - ou + ?
                else              bi -= ds*Y[i+1] + di*Y[i-1];

                Y[i] = bi/d; //calcula incognita

            }
        
        
        }
    
    }

}




int main(){

    // Edo *edoeq = (Edo *)malloc(sizeof(Edo));
    // //y'' = 6x - 0.5x², x e (0,12)
    // //y(0)=0 y(12) = 0
    // edoeq->a = 0;
    // edoeq->b = 12;
    // edoeq->n = 5;
    // edoeq->ya = 0;
    // edoeq->yb = 0;
    // edoeq->p = getP;
    // edoeq->q = getQ;
    // edoeq->r = getR;
    // double *Y;
    // Y = (double *)calloc((edoeq->n), sizeof(double));
    // gaussSeidel(edoeq, Y);
    // for(int i=0; i<edoeq->n; i++)
    //     printf("%f ", Y[i]);
    
    //****************************************************
    Edo2 *edoeq = (Edo2 *)malloc(sizeof(Edo2));
    edoeq->Lx = M_PI;
    edoeq->Ly = M_PI_2;
    edoeq->n = 5;
    edoeq->m = 3;
    edoeq->func = func;
    edoeq->u1 = u1;
    edoeq->u2 = u2;
    edoeq->u3 = u3;
    edoeq->u4 = u4;



}