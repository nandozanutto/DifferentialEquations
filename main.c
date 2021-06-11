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
    double (* u1)(double), (* u2)(double), (* u3)(double), (* u4)(double); //(0,y) (lx, y) (x, 0) (x, ly)
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
    //sin²(x)
    return sin(x)*sin(x);
}

double u1(double var){
    return 20.0;
}

double u2(double var){
    return 45.0;
}

double u3(double var){
    return 0;
}

double u4(double var){
    return 100;
}


double gaussSeidel(Edo *edoeq, double *Y)
{
    int n = edoeq->n, k, i;
    double h, xi, bi, yj, d, di, ds;
    double r[n], norma=0;
    double biV[n], diV[n], dV[n], dsV[n];

    h = (edoeq->b - edoeq->a)/(n+1); //Largura do passo da malha
    for(k=0; k<50; ++k){
        for(i=0; i<n; ++i){
            xi = edoeq->a + (i+1)*h;//valor xi da malha
            bi = h*h * edoeq->r(xi);//termo independente
            biV[i] = bi;
            di = 1 - h*edoeq->p(xi)/2.0;//diagonal inferior
            d = -2 + h*h * edoeq->q(xi);//diagonal principal
            ds = 1 + h*edoeq->p(xi)/2.0;//diagonal superior
            if(i==0)          bi -= ds*Y[i+1] + edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0);
            else if(i == n-1) bi -= di*Y[i-1] + edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0);
            else              bi -= ds*Y[i+1] + di*Y[i-1];

            Y[i] = bi/d; //calcula incognita
            diV[i] = di;
            dV[i] = d;
            dsV[i] = ds;
        }

    }

    for(int i=0; i<n; i++){
        if(i==0) r[i] = biV[i] - (dsV[i]*Y[i+1] + dV[i]*Y[i] + edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0));
        else if(i == n-1)  r[i] = biV[i] - (diV[i]*Y[i-1] + dV[i]*Y[i] + edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0));
        else r[i] = biV[i] - (diV[i]*Y[i-1] + dV[i]*Y[i] + dsV[i]*Y[i+1]);
    }

    for(int i=0; i<n; i++)
        norma += r[i]*r[i];
    norma = sqrt(norma);
    return norma;

}

int gaussSeidel2(Edo2 *edoeq)
{
    int n = edoeq->n, m = edoeq->m, i, k, j;
    double U[n][m], r[n][m], norma=0;
    for(int j=0; j<edoeq->m; j++)//zerando a matriz
        for(int i=0; i<edoeq->n; i++)
            U[i][j] = 0;


    double hx, hy, xi, bi, yj, d, di, ds, di2, ds2, biV[n][m];

    hx = edoeq->Lx/(n+1);
    hy = edoeq->Ly/(m+1);//largura do passo
    
    for(k=0; k<50; ++k){
        
        for(j=0; j<m; j++){
            
            for(i=0; i<n; i++){
                xi = (i+1)*hx;//valor xi da malha
                yj = (j+1)*hy;//valor yj da malha

                bi = hx*hx*hy*hy*edoeq->func(xi, yj);//termo indep.
                biV[i][j] = bi;
                di = hy*hy;
                di2 = hx*hx;
                d = -2.0*(hx*hx + hy*hy) - 1;//?????
                ds = hy*hy;
                ds2 = hx*hx;


                if((i == 0) && (j == 0)) bi -= di2*edoeq->u3(xi) + di*edoeq->u1(yj) + ds*U[i+1][j] + ds2*U[i][j+1];
                else if((i == 0) && (j != m-1)) bi -= di2*U[i][j-1] + di*edoeq->u1(yj) + ds*U[i+1][j] + ds2*U[i][j+1];
                else if((i == 0) && (j == m-1)) bi -= di2*U[i][j-1] + di*edoeq->u1(yj) + ds*U[i+1][j] + ds2*edoeq->u4(xi);
                
                else if((i == n-1) && (j == 0)) bi -= di2*edoeq->u3(xi) + di*U[i-1][j] + ds*edoeq->u2(yj) + ds2*U[i][j+1];
                else if((i == n-1) && (j != m-1)) bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*edoeq->u2(yj) + ds2*U[i][j+1];
                else if((i == n-1) && (j == m-1)) bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*edoeq->u2(yj) + ds2*edoeq->u4(xi);
                
                else if((i != n-1) && (j == 0)) bi -= di2*edoeq->u3(xi) + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1];
                else if((i != n-1) && (j != m-1)) bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1];
                else if((i != n-1) && (j == m-1))  bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] +  ds2*edoeq->u4(xi); 
                U[i][j] = bi/d; //calcula incognita

            }

        }
        

    
    }

    for(j=0; j<m; j++){
        for(i=0; i<n; i++){
            xi = (i+1)*hx;//valor xi da malha
            yj = (j+1)*hy;//valor yj da malha
            if((i == 0) && (j == 0)) r[i][j] = biV[i][j] - (di2*edoeq->u3(xi) + di*edoeq->u1(yj) + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            else if((i == 0) && (j != m-1)) r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*edoeq->u1(yj) + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            else if((i == 0) && (j == m-1)) r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*edoeq->u1(yj) + ds*U[i+1][j] + ds2*edoeq->u4(xi) + d*U[i][j]);
            
            else if((i == n-1) && (j == 0)) r[i][j] = biV[i][j] - (di2*edoeq->u3(xi) + di*U[i-1][j] + ds*edoeq->u2(yj) + ds2*U[i][j+1] + d*U[i][j]);
            else if((i == n-1) && (j != m-1)) r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*edoeq->u2(yj) + ds2*U[i][j+1] + d*U[i][j]);
            else if((i == n-1) && (j == m-1)) r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*edoeq->u2(yj) + ds2*edoeq->u4(xi) + d*U[i][j]);
            
            else if((i != n-1) && (j == 0)) r[i][j] = biV[i][j] - (di2*edoeq->u3(xi) + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            else if((i != n-1) && (j != m-1)) r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            else if((i != n-1) && (j == m-1))  r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] +  ds2*edoeq->u4(xi) + d*U[i][j]);
            printf("r: %f ", r[i][j]);
        }
    }
    for(j=0; j<m; j++)
        for(i=0; i<n; i++)
            norma+=r[i][j]*r[i][j];
    norma = sqrt(norma);
    for(int j=0; j<edoeq->m; j++){
        for(int i=0; i<edoeq->n; i++)
            printf("%.7g ", U[i][j]);
        printf("\n");
    }
    printf("%1.27g", norma);
}


int main(){

    // Edo *edoeq = (Edo *)malloc(sizeof(Edo));
    // //y'' = 6x - 0.5x², x e (0,12)
    // //y(0)=0 y(12) = 0
    // edoeq->a = 0;
    // edoeq->b = 12;
    // edoeq->n = 10;
    // edoeq->ya = 0;
    // edoeq->yb = 0;
    // edoeq->p = getP;
    // edoeq->q = getQ;
    // edoeq->r = getR;
    // double *Y;
    // Y = (double *)calloc((edoeq->n), sizeof(double));
    // gaussSeidel(edoeq, Y);
    // for(int i=0; i<edoeq->n; i++)
    //     printf("%.7g ", Y[i]);
    
    //****************************************************
    Edo2 *edoeq = (Edo2 *)malloc(sizeof(Edo2));
    edoeq->Lx = 6;
    edoeq->Ly = 8;
    edoeq->n = 5;
    edoeq->m = 3;
    edoeq->func = func;
    edoeq->u1 = u1;
    edoeq->u2 = u2;
    edoeq->u3 = u3;
    edoeq->u4 = u4; 
    printf("%d",gaussSeidel2(edoeq));


}