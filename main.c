#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

/*  Retorna tempo em milisegundos
    Forma de uso:
 
    real_t tempo;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/

typedef float real_t;

real_t timestamp(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((real_t)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}


//Matriz tridiagonal
typedef struct {
    int n; //numero de pontos internos na malha
    real_t a, b; //intervalo
    real_t ya, yb;//condições de contorno
    real_t (* p)(real_t), (* q)(real_t, int), (* r)(real_t, int);
}Edo;

typedef struct {
    int n, m; //numero de pontos internos na malha
    real_t Lx, Ly; //intervalo
    real_t (* u1)(real_t, int), (* u2)(real_t, int), (* u3)(real_t, int), (* u4)(real_t, int); //(0,y) (lx, y) (x, 0) (x, ly)
    real_t (* func)(real_t, real_t, int);
}Edo2;

real_t getP(real_t x){
    return 0;
}

real_t getQ(real_t x, int isFirst){
    if(isFirst) return 0;
    else return 1;
}

real_t getR(real_t x, int isFirst){//Eq a.
    //6x - 0.5x²
    if(isFirst) return 6*x - 0.5*x*x;
    else return 0;
}

real_t func(real_t x, real_t y, int isFirst){//Eq a.
    //sin²(x)
    if(isFirst) return sin(x)*sin(x);
    else return -cos(x+y) -cos(x-y);
}

real_t u1(real_t var, int isFirst){
    if(isFirst) return 20.0;
    else return cos(var);
}

real_t u2(real_t var, int isFirst){
    if(isFirst) return 45.0;
    else return -cos(var);
}

real_t u3(real_t var, int isFirst){
    if(isFirst) return 0;
    else return cos(var);
}

real_t u4(real_t var, int isFirst){
    if(isFirst) return 100;
    else return 0;
}

void printV(real_t *v, int n){
    for(int i=0; i<n; i++)
        printf("%.7g ", v[i]);
    printf("\n");
}

real_t gaussSeidel(Edo *edoeq, real_t *Y, int isFirst)
{
    real_t tTotal = timestamp();
    printf("%f", tTotal);
    int n = edoeq->n, k, i;
    real_t h, xi, bi, yj, d, di, ds;
    real_t r[n], norma=0;
    real_t biV[n], diV[n], dV[n], dsV[n];

    h = (edoeq->b - edoeq->a)/(n+1); //Largura do passo da malha
    for(k=0; k<50; ++k){
        for(i=0; i<n; ++i){
            xi = edoeq->a + (i+1)*h;//valor xi da malha
            bi = h*h * edoeq->r(xi, isFirst);//termo independente
            biV[i] = bi;
            di = 1 - h*edoeq->p(xi)/2.0;//diagonal inferior
            d = -2 + h*h * edoeq->q(xi, isFirst);//diagonal principal
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
    tTotal = timestamp() - tTotal;
    for(int i=0; i<n; i++){//fazendo residuo
        if(i==0) r[i] = biV[i] - (dsV[i]*Y[i+1] + dV[i]*Y[i] + edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0));
        else if(i == n-1)  r[i] = biV[i] - (diV[i]*Y[i-1] + dV[i]*Y[i] + edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0));
        else r[i] = biV[i] - (diV[i]*Y[i-1] + dV[i]*Y[i] + dsV[i]*Y[i+1]);
    }

    for(int i=0; i<n; i++)//fazendo norma
        norma += r[i]*r[i];
    norma = sqrt(norma);


    if(isFirst) printf("***** item (a): n = %d, H = %.7g\nSL:\n", n, h);
    else printf("***** item (b): n = %d, H = %.7g\nSL:\n", n, h);
    printV(dsV, n);
    printV(dV, n);
    printV(diV, n);
    printV(biV, n);
    printf("Y: ");
    printV(Y, n);
    printf("Norma L2: %.7g, Tempo: %f ms,\n", norma, tTotal);

    return norma;

}

int gaussSeidel2(Edo2 *edoeq, int isFirst)
{
    real_t tTotal = timestamp();
    int n = edoeq->n, m = edoeq->m, i, k, j;
    real_t U[n][m], r[n][m], norma=0;
    for(int j=0; j<edoeq->m; j++)//zerando a matriz
        for(int i=0; i<edoeq->n; i++)
            U[i][j] = 0;


    real_t hx, hy, xi, bi, yj, d, di, ds, di2, ds2, biV[n][m];

    hx = edoeq->Lx/(n+1);
    hy = edoeq->Ly/(m+1);//largura do passo
    
    for(k=0; k<50; ++k){
        
        for(j=0; j<m; j++){
            
            for(i=0; i<n; i++){
                xi = (i+1)*hx;//valor xi da malha
                yj = (j+1)*hy;//valor yj da malha

                bi = hx*hx*hy*hy*edoeq->func(xi, yj, isFirst);//termo indep.
                biV[i][j] = bi;
                di = hy*hy;
                di2 = hx*hx;
                d = isFirst ? -2.0*(hx*hx + hy*hy) - 1 : -2.0*(hx*hx + hy*hy);
                ds = hy*hy;
                ds2 = hx*hx;


                if((i == 0) && (j == 0)) bi -= di2*edoeq->u3(xi, isFirst) + di*edoeq->u1(yj, isFirst) + ds*U[i+1][j] + ds2*U[i][j+1];
                else if((i == 0) && (j != m-1)) bi -= di2*U[i][j-1] + di*edoeq->u1(yj, isFirst) + ds*U[i+1][j] + ds2*U[i][j+1];
                else if((i == 0) && (j == m-1)) bi -= di2*U[i][j-1] + di*edoeq->u1(yj, isFirst) + ds*U[i+1][j] + ds2*edoeq->u4(xi, isFirst);
                
                else if((i == n-1) && (j == 0)) bi -= di2*edoeq->u3(xi, isFirst) + di*U[i-1][j] + ds*edoeq->u2(yj, isFirst) + ds2*U[i][j+1];
                else if((i == n-1) && (j != m-1)) bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*edoeq->u2(yj, isFirst) + ds2*U[i][j+1];
                else if((i == n-1) && (j == m-1)) bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*edoeq->u2(yj, isFirst) + ds2*edoeq->u4(xi, isFirst);
                
                else if((i != n-1) && (j == 0)) bi -= di2*edoeq->u3(xi, isFirst) + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1];
                else if((i != n-1) && (j != m-1)) bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1];
                else if((i != n-1) && (j == m-1))  bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] +  ds2*edoeq->u4(xi, isFirst); 
                U[i][j] = bi/d; //calcula incognita

            }

        }
        

    
    }

    for(j=0; j<m; j++){//fazendo residuo
        for(i=0; i<n; i++){
            xi = (i+1)*hx;//valor xi da malha
            yj = (j+1)*hy;//valor yj da malha
            if((i == 0) && (j == 0)) r[i][j] = biV[i][j] - (di2*edoeq->u3(xi, isFirst) + di*edoeq->u1(yj, isFirst) + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            else if((i == 0) && (j != m-1)) r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*edoeq->u1(yj, isFirst) + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            else if((i == 0) && (j == m-1)) r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*edoeq->u1(yj, isFirst) + ds*U[i+1][j] + ds2*edoeq->u4(xi, isFirst) + d*U[i][j]);
            
            else if((i == n-1) && (j == 0)) r[i][j] = biV[i][j] - (di2*edoeq->u3(xi, isFirst) + di*U[i-1][j] + ds*edoeq->u2(yj, isFirst) + ds2*U[i][j+1] + d*U[i][j]);
            else if((i == n-1) && (j != m-1)) r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*edoeq->u2(yj, isFirst) + ds2*U[i][j+1] + d*U[i][j]);
            else if((i == n-1) && (j == m-1)) r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*edoeq->u2(yj, isFirst) + ds2*edoeq->u4(xi, isFirst) + d*U[i][j]);
            
            else if((i != n-1) && (j == 0)) r[i][j] = biV[i][j] - (di2*edoeq->u3(xi, isFirst) + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            else if((i != n-1) && (j != m-1)) r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            else if((i != n-1) && (j == m-1))  r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] +  ds2*edoeq->u4(xi, isFirst) + d*U[i][j]);
        }
    }
    for(j=0; j<m; j++) //fazendo norma
        for(i=0; i<n; i++)
            norma+=r[i][j]*r[i][j];
    norma = sqrt(norma);
    for(int j=0; j<edoeq->m; j++){
        for(int i=0; i<edoeq->n; i++)
            printf("%.7g ", U[i][j]);
        printf("\n");
    }
    printf("norma = %.7g", norma);
    tTotal = timestamp() - tTotal;
}


int main(){

    Edo *edoeq = (Edo *)malloc(sizeof(Edo));
    //y'' = 6x - 0.5x², x e (0,12)
    //y(0)=0 y(12) = 0
    edoeq->a = 0;
    edoeq->b = 12;
    edoeq->n = 10;
    edoeq->ya = 0;
    edoeq->yb = 0;
    edoeq->p = getP;
    edoeq->q = getQ;
    edoeq->r = getR;
    real_t *Y;
    Y = (real_t *)calloc((edoeq->n), sizeof(real_t));
    gaussSeidel(edoeq, Y, 1);
    
    printf("\n******************************************\n");
    //************************************************
    
    // Edo *edoeq = (Edo *)malloc(sizeof(Edo));
    // //y'' + y = 0
    // //y(0)=0 y(12) = 0
    // edoeq->a = 0;
    // edoeq->b = 1;
    // edoeq->n = 10;
    // edoeq->ya = 0;
    // edoeq->yb = 1;
    // edoeq->p = getP;
    // edoeq->q = getQ;
    // edoeq->r = getR;
    // real_t *Y;
    // Y = (real_t *)calloc((edoeq->n), sizeof(real_t));
    // gaussSeidel(edoeq, Y, 0);
    // for(int i=0; i<edoeq->n; i++)
    //     printf("%.7g ", Y[i]);
    
    //****************************************************
    Edo2 *edoeq2 = (Edo2 *)malloc(sizeof(Edo2));
    edoeq2->Lx = 6;
    edoeq2->Ly = 8;
    edoeq2->n = 5;
    edoeq2->m = 3;
    edoeq2->func = func;
    edoeq2->u1 = u1;
    edoeq2->u2 = u2;
    edoeq2->u3 = u3;
    edoeq2->u4 = u4; 
    printf("%d",gaussSeidel2(edoeq2, 1));

    // ****************************************************
    // Edo2 *edoeq = (Edo2 *)malloc(sizeof(Edo2));
    // edoeq->Lx = M_PI;
    // edoeq->Ly = M_PI_2;
    // edoeq->n = 5;
    // edoeq->m = 3;
    // edoeq->func = func;
    // edoeq->u1 = u1;
    // edoeq->u2 = u2;
    // edoeq->u3 = u3;
    // edoeq->u4 = u4; 
    // printf("%d",gaussSeidel2(edoeq, 0));


}