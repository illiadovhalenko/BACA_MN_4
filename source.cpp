//Illia Dovhalenko
#include "vectalg.h"
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;


Vector solveEquations(
        const Matrix & A,   // Macierz
        const Vector & b,   // Wektor
        double  eps         // dopuszczalny błąd
){
    Matrix Aa=A;
    Vector bb=b;
    Vector odpowiedz(b.size());
    Vector normy_wierszow(b.size());
    vector<int> porzadek(b.size());
    for(int i=0; i< porzadek.size(); i++){
        porzadek[i]=i;
    }
    for(int k=0; k<Aa.size()-1; k++){
        double max=0;
        int max_j=0;
        for(int j=k; j<Aa.size(); j++){
            if(abs(Aa(porzadek[j], k))>max){
                max=abs(Aa(porzadek[j], k))/normy_wierszow[j];
                max_j=j;
            }
        }
        swap(porzadek[k], porzadek[max_j]);
        for (int i = k+1; i < Aa.size(); i++) {
            double temp = Aa(porzadek[i], k) / Aa(porzadek[k], k);
            Aa(porzadek[i], k)=temp;
            for (int j = k + 1; j < Aa.size(); j++)
                Aa(porzadek[i], j) -= Aa(porzadek[k], j) * temp;
            bb[porzadek[i]]-=Aa(porzadek[i], k)*bb[porzadek[k]];
        }
    }
    for(int i=A.size()-1; i>=0; i--){
        double suma=0;
        for(int j=i+1; j<A.size(); j++)
            suma+=Aa(porzadek[i], j)*odpowiedz[j];
        odpowiedz[i]=(bb[porzadek[i]]-suma)/Aa(porzadek[i], i);
    }
    Vector rk(odpowiedz.size());
    do{
        rk = residual_vector(A, b, odpowiedz);
        if(rk.max_norm()<eps)
            return odpowiedz;
        transform(odpowiedz.begin(), odpowiedz.end(), solveEquations(A, rk, eps).begin(), odpowiedz.begin(), plus<long double>());
    }while(true);
}