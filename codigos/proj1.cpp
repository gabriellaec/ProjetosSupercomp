#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include<bits/stdc++.h>

using namespace std;
using std::max;


struct item {
    int valor;
    int salto_i;
    int salto_j;
};


int w(string a, string b);
item find_max(int diagonal, int delecao, int insercao);
int smith_waterman(vector<string> a, vector<string> b, int n, int m);


int main() {
    int m=0;
    int n=0;
    vector<string> a;
    vector<string> b;
    int max_score;

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

// Imprimido sequências iniciais
    cout << endl << "----- Sequências iniciais -----" << endl;
    cout << "a - tamanho: " << n << endl;
    for(auto& el: a){
        cout << el;
    }
    cout << endl;
    cout << "b - tamanho: " << m << endl;
    for(auto& el: b){
        cout << el;
    }

// Obtenção do score e das sequências alinhadas
    max_score = smith_waterman(a,b,n,m);

// Score obtido
    cout << endl << endl << "----- Score -----" << endl;
    cout << "max_score: "<< max_score << endl << endl;

     return 0;
}


int w(string a, string b){
    if (a == b) return 2;  //match
    else return -1;  //mismatch or gap
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


int smith_waterman(vector<string> a, vector<string> b, int n, int m){
    item H[n+1][m+1];
    

    // Zerando colunas especificadas
    for (int i = 0; i <= n; i++) {
        H[i][0].valor=0;
    }
    for (int j = 0; j <= m; j++) {
        H[0][j].valor = 0;
    }


    // Obtendo matriz e achando o valor máximo
    int diagonal, delecao, insercao;
    int maximo_H = 0; int max_val_i = 0; int max_val_j = 0;

    for (int i=1; i<=n; i++){
        for (int j=1; j<=m; j++){
            diagonal = H[i-1][j-1].valor + w(a[i-1],b[j-1]);
            delecao = H[i-1][j].valor - 1;
            insercao = H[i][j-1].valor - 1;

            H[i][j]=find_max(diagonal, delecao, insercao);

            if (H[i][j].valor > maximo_H) {
                    maximo_H=H[i][j].valor;
                    max_val_i=i;
                    max_val_j=j;
            }
        }
    }


    // Imprimindo a matriz obtida
    cout << endl << endl << "----- Matriz -----" << endl;
    for (int i=0; i<=n; i++){
        for (int j=0; j<=m; j++){
            cout << H[i][j].valor << " ";
        }
        cout << endl;
    }


// ____________________________________________ //
//  Reconstrução do alinhamento das sequências //
// __________________________________________ //

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

    // Imprimindo resultado
    cout << endl << "----- Reconstrução do alinhamento -----" << endl;
    cout << "a: ";
    for(auto& el: match_seq_a){
        cout << el;
    }
    cout << endl;
    cout << "b: ";
    for(auto& el: match_seq_b){
        cout << el;
    }


    return maximo_H;
}

// Para compilar: 

// g++ -Wall -O3 proj1.cpp -o proj1 
// ./proj1 < dna.seq