#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include<bits/stdc++.h>

using namespace std;
using std::max;

#define WMAT 2
#define WMIS -1
#define WGAP -1

struct item {
    int valor;
    int salto_i;
    int salto_j;
};

struct resultado {
    int valor;
    vector<string> seq1;
    vector<string> seq2;
};


int w(string a, string b);
item find_max(int diagonal, int delecao, int insercao);
resultado smith_waterman_results(vector<string> a, vector<string> b, int n, int m);


int main() {
    int m=0;
    int n=0;
    vector<string> a;
    vector<string> b;

 // Lendo tamanho das sequências   
    cin >> n >> m;

// Lendo sequências iniciais
    string base;
    cin >> base;
    for(int i = 0; i < n; i++){
        a.push_back({base[i]});
    }
    cin >> base;
    for(int i = 0; i < m; i++){
        b.push_back({base[i]});
    }

    resultado result =smith_waterman_results(a,b,n,m);

    cout << "--- valor maximo ---" << endl;
    cout << result.valor << endl;
    
    cout << endl << "----- Reconstrução do alinhamento -----" << endl;
    cout << "a: ";
    for(auto& el: result.seq1){
        cout << el;
    }
    cout << endl;
    cout << "b: ";
    for(auto& el: result.seq2){
        cout << el;
    }
    cout << endl << endl;

     return 0;
}


int w(string a, string b){
    if (a == b)
        return WMAT;  //match
    else if (a != b)
        return WMIS;
   else
        return WGAP;  // gap
}


item find_max(int diagonal, int delecao, int insercao){
    item maior;
    int max_score_local=0;
    max_score_local = max({0, diagonal, delecao, insercao});

    maior.valor = 0;
    maior.salto_i = 0;
    maior.salto_j = 0;
    if (diagonal == max_score_local){
        maior.valor = diagonal;
        maior.salto_i = 1;  // salto na diagonal
        maior.salto_j = 1;
    } 

    else if (delecao == max_score_local){
        maior.valor = delecao;
        maior.salto_i = 1;  // salto para cima
        maior.salto_j = 0;
    } 

    else if (insercao == max_score_local){
        maior.valor = insercao;
        maior.salto_i = 0;
        maior.salto_j = 1; // salto da esquerda para a direita
    } 

    return maior;
}



resultado smith_waterman_results(vector<string> a, vector<string> b, int n, int m){
    vector<vector<item>> H;
    int maximo_H = 0; int max_val_i = 0; int max_val_j = 0;

    // Obtendo matriz e achando o valor máximo
    int diagonal, delecao, insercao;

    H.resize(n+1);
    for (int i=0; i<=n; i++){
        H[i].resize(m+1);
    }
    for (int i=1; i<=n; i++){
        for (int j=1; j<=m; j++){
            diagonal = H[i-1][j-1].valor + w(a[i-1],b[j-1]);
            delecao = H[i-1][j].valor - 1;
            insercao = H[i][j-1].valor - 1;

            H[i][j]=(item)(find_max(diagonal, delecao, insercao));

            if (H[i][j].valor > maximo_H) {
                    maximo_H=H[i][j].valor;
                    max_val_i=i;
                    max_val_j=j;
            }

        }
    }



    vector<string> match_seq_a;
    vector<string> match_seq_b;
    int i=max_val_i; int j=max_val_j;

    // Reconstruindo o caminho a partir dos saltos do struct
    while ( (i>0 && j>0)  && (!(H[i][j].salto_j==0 && H[i][j].salto_i==0)) ) {
        int pos_i=i;
        int pos_j=j;
        if (H[i][j].valor == 0) break; // célula da matriz com valor zero

        if (H[pos_i][pos_j].salto_i==0 && H[pos_i][pos_j].salto_j ==1){
            match_seq_a.push_back("_");
            match_seq_b.push_back(b[j-1]);
        }
        else if (H[pos_i][pos_j].salto_i==1 && H[pos_i][pos_j].salto_j ==0){
            match_seq_a.push_back(a[i-1]);
            match_seq_b.push_back("_");
        }
        else{
            match_seq_a.push_back(a[i-1]);
            match_seq_b.push_back(b[j-1]);
        }
        
        
        i= i- H[pos_i][pos_j].salto_i;
        j=j- H[pos_i][pos_j].salto_j;                
    }

    // Invertendo sequências
    reverse(match_seq_a.begin(),match_seq_a.end());
    reverse(match_seq_b.begin(),match_seq_b.end());
    
    return {maximo_H, match_seq_a, match_seq_b};
}


// Para compilar: 

// g++ -Wall -O3 heuristica.cpp -o heuristica 
// ./heuristica < dna.seq > out.txt