#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include<bits/stdc++.h>

using namespace std;
using std::max;


// ----- Pesos definidos ----- //
#define WMAT 2
#define WMIS -1
#define WGAP -1


// ----- Structs ----- //
struct item {
    int item_score;
    vector<string> seq_a;
    vector<string> seq_b;
};

struct item_sw{
    int valor;
    int salto_i;
    int salto_j;
};


// ----- Declarações de funções ----- //
int m,n;
int w(string a, string b);
void gera_subseq(vector<string> seq, int start_point, int end_point, vector<vector<string>>& matriz_subseq);
int score_simples(vector<string> sa, vector<string> sb);



int main() {
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

    vector<vector<string>> subseqs_a;
    vector<vector<string>> subseqs_b;
    gera_subseq(a, 0, 1, subseqs_a);
    gera_subseq(b, 0, 1, subseqs_b);

// ----- Obtendo o score e as sequências alinhadas ----- //
    vector<string> melhor_sub_a, melhor_sub_b, subvector_a,subvector_b;
    int score=0, max_score=0;

    for (auto& sub_a : subseqs_a){
        for (auto& sub_b : subseqs_b){
            vector<string> sub_a_atual, sub_b_atual;
            if ((int)sub_a.size() == (int)sub_b.size()){
                score = score_simples(sub_a, sub_b);
                sub_a_atual=sub_a;
                sub_b_atual=sub_b;
            }
            else if ((int)sub_a.size() > (int)sub_b.size()){
                vector<string> subvector_a = {sub_a.begin(), sub_a.begin()+(int)sub_b.size()}; 
                score = score_simples(subvector_a, sub_b);
                sub_a_atual=subvector_a;
                sub_b_atual=sub_b;
                
            }else{
                vector<string> subvector_b = {sub_b.begin(), sub_b.begin()+(int)sub_a.size()}; 
                score = score_simples(sub_a, subvector_b);
                sub_a_atual=sub_a;
                sub_b_atual=subvector_b;
            }
        
            if (score>max_score){
                max_score=score;
                melhor_sub_a = sub_a_atual;
                melhor_sub_b = sub_b_atual;
            }
        }
    }
    

// ----- Imprimindo Score obtido ----- //
    cout << "----- Score -----" << endl;
    cout << "max_score: "<< max_score << endl << endl;

    cout << " ----- Melhor subsequência de A -----" << endl;
    for (auto& el : melhor_sub_a ){
        cout << el << " ";
    }
    cout << endl;

    cout << " ----- Subsequência B correspondente -----" << endl;
    for (auto& el : melhor_sub_b){
        cout << el << " ";
    }
    cout << endl;
    

     return 0;
}


// ----- Função que gera todas as subsequências possíveis ----- //
void gera_subseq(vector<string> seq, int start_point, int end_point, vector<vector<string>>& matriz_subseq){
    if (end_point > (int)seq.size())
      return;
    else if (start_point > end_point){
        gera_subseq(seq, 0, end_point+1, matriz_subseq);
    }else{
        if (start_point != end_point){
            vector<string> subseq;
            for (int j=start_point; j<end_point; j++){
                subseq.push_back(seq[j]);
            }
            matriz_subseq.push_back(subseq);
        }
        gera_subseq(seq, start_point+1, end_point, matriz_subseq);
    }
}


// ----- Função que calcula o score do alinhamento simples entre duas sequências ----- //
int score_simples(vector<string> sa, vector<string> sb){
    int score=0; 
    vector<string> melhor_subseq_a, melhor_subseq_b;
    for (int i=0; i<(int)sa.size(); i++){
        score += w(sa[i],sb[i]); 
    }  

    return score;
}


// ----- Função que retorna os pesos do cálculo do score ----- //
int w(string a, string b){
    if (a == b)
        return WMAT;  //match
    else if (a != b)
        return WMIS;  // mismatch
   else
        return WGAP;  // gap

}



// Para compilar: 

// g++ -Wall -O3 p3_simple_score.cpp -o p3_simple_score
// ./p3_simple_score < dna.seq