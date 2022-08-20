// Written by Robert Gerbicz

// multi threaded version, use the the -t switch, say -t 8 to run it on 8 threads, if not set, then using all threads.

// my long compilation line:
// g++ -mtune=skylake -march=skylake -fomit-frame-pointer -flto -frename-registers -mavx2 -m64 -O2 -o s 5x5multi.c -fopenmp

// Fast search for Matt Parker's 5 five letter words that has 25 different letters from a-z.

// in nutshell the algorithm:
// define a new alphabet where rare letters comes first, sort the words by this way.
// And use that there is only one missing letter, it means that
// at each depth there could be at most 2 letters what the word should use in the dfs algorithm.
// Say the first word should use q or j as a letter. And so on, but this 1/2 
// letters could depend on the previous words (since they can contain any letter).
// Furthermore on each depth compute what words could be on the remaining depths,
// for this use the results from the previous depth (memo).

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <sys/time.h>
#include <cstring>
#include "omp.h"// for multi thread

using namespace std;

#define NL 5 // number of letters on each word
#define ASIZE 26 // size of alphabet

int numwords;
vector<pair<int,int>>w;
unordered_map<int,vector<string>>mymap;

void inits(void){
    
    assert(ASIZE%NL==1);

    ifstream fin("vocabulary.txt");
    if(!fin){
        cout<<"The words should be in vocabulary.txt file (one word per line using 'a'-'z' letters)"<<endl;
        exit(1);
    }

    int i,j,k;
    string st;
    vector<string>u;
    unordered_set<int>myset;

    w.clear();
    mymap.clear();
    
    vector<pair<int,int>>freq(26,{0,0});    
    for(i=0;i<26;i++)freq[i].second=i;
    
    numwords=0;
    while(getline(fin,st)){
        st.erase(std::remove_if(st.begin(),st.end(),[](char c){return !std::isalnum(c); }),st.end());// to remove non-letters
        
        if(st.length()!=NL) continue;
        
        numwords++;
        
        int hash=0;
        bool cont=false;
        for(i=0;i<NL;i++){
            int pos=st[i]-'a';
            if((hash>>pos)&1){cont=true;break;}
            hash+=(1<<pos);
        }
        if(!cont){
            u.push_back(st);
            if(myset.count(hash)==0){
               for(i=0;i<NL;i++)freq[st[i]-'a'].first++;
               myset.emplace(hash);
            }
        }        
    }
    fin.close();
    
    sort(freq.begin(),freq.end());
    
    int inv[ASIZE];
    for(i=0;i<ASIZE;i++)inv[freq[i].second]=i;
  
    myset.clear();
    
    for(i=0;i<u.size();i++){
        int h0=0,h1=0;
        st=u[i];
        for(j=0;j<NL;j++){
            k=st[j]-'a';
            h0+=(1<<(ASIZE-1-inv[k]));
            h1+=(1<<inv[k]);
        }
        
        mymap[h1].push_back(u[i]);
        if(myset.count(h1)==0){
            myset.emplace(h1);
            w.push_back({h0,h1});
        }
    }
    sort(w.begin(),w.end());
    reverse(w.begin(),w.end());
}


void fun(int nthreads,vector<vector<string>>&sols){
    
    sols.clear();
    int i;
    int th_id=omp_get_thread_num();
    int maxdepth=(ASIZE-1)/NL;
    int nw=w.size();
    int T[nw*(maxdepth+1)];
    for(i=0;i<nw;i++)
        T[i]=w[i].second;
    
    int depth=1;
    int excess[NL+1];// on each depth: excess=0,1 
    int P[NL+1],offset[NL+1],shift[NL+1];
    P[1]=-1;
    excess[0]=0;// excess[depth-1]=1 iff we already skipped a letter, otherwise it is 0
    offset[0]=0;// offset[depth-1] to store the used letters bit offsets
    shift[0]=0; // shift[depth-1] gives the first letter that we should consider
    
    int ntrace=0;
    int num_solution=0;
    int startpos[NL+1],endpos[NL+1];
    startpos[1]=0;
    endpos[1]=nw;
    
    while(depth>0){
        ntrace++;
        P[depth]++;
        
        if(P[depth]>=endpos[depth]){
            depth--;
            continue;
        }
        
        if(depth==1&&P[depth]%nthreads!=th_id)  continue;        
        
        if(depth==maxdepth){
           // solution found!!! save, the anagrams also
           int depth2=1;
           int P2[maxdepth+1];
           int sizes[maxdepth+1];
           for(i=1;i<=maxdepth;i++)sizes[i]=mymap[T[P[i]]].size();
           P2[1]=-1;
           while(depth2>0){
               P2[depth2]++;
               if(P2[depth2]>=sizes[depth2]){depth2--;continue;}
               if(depth2==maxdepth){
                   vector<string>a;
                   for(i=1;i<=maxdepth;i++)
                       a.push_back(mymap[T[P[i]]][P2[i]]);
                   sols.push_back(a);
                   continue;
               }
               depth2++;
               P2[depth2]=-1;
           }
           continue;   
        }
        
        int hv=(offset[depth-1]|T[P[depth]]);
        offset[depth]=hv;
        
        int tmp=(hv>>shift[depth-1]);
        excess[depth]=excess[depth-1];
        int numones=0;
        while(tmp&1){
            tmp>>=1;
            numones++;
        }
        
        if(numones==0){
            tmp>>=1;
            excess[depth]++;
            numones++;
            if((tmp&1)==0||excess[depth]>1){
                // already at least two characters are missing, backtrack!!!
                depth--;
                continue;
            }
            while(tmp&1){
                tmp>>=1;
                numones++;
            }
        }
        shift[depth]=shift[depth-1]+numones;
        
        int size2=0,position=endpos[depth];
        for(i=P[depth]+1;i<endpos[depth];i++){
            if((hv & T[i])==0){
                T[position++]=T[i];
                size2++;
            }
        }

        if(size2+depth>=maxdepth){
            depth++;
            startpos[depth]=endpos[depth-1];
            endpos[depth]=position;
            P[depth]=startpos[depth]-1;
        }
    }
}

void solve(int nthreads){
    
    struct timeval start, end;  
    gettimeofday(&start, NULL);
    double time_taken;

    int maxthreads=omp_get_max_threads();   
    if(nthreads==0)nthreads=maxthreads;
    assert(nthreads>=0&&nthreads<=maxthreads);   
    omp_set_num_threads(nthreads);    

    inits();//it is single threaded
    
    vector<vector<vector<string>>>ans(nthreads);
    
    #pragma omp parallel for schedule(dynamic,1)
    for(int i=0;i<nthreads;i++)
        fun(nthreads,ans[i]);
    
    int count=0;
    for(int h=0;h<nthreads;h++){
    for(int i=0;i<ans[h].size();i++){
        cout<<++count<<" ";
        for(auto a:ans[h][i])cout<<a<<" ";
        cout<<endl;
    }}
    
    gettimeofday(&end, NULL); 
    time_taken = (end.tv_sec - start.tv_sec) * 1e6;
    time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
    printf("Computed all solutions on %d threads, Wall time taken in seconds = %lf.\n",nthreads,time_taken);
    printf("Number of words=%d, the pruned size=%d.\n",numwords,(int)w.size());
}

int main(int argc, char *argv[]){
    
   int nthreads=0;
   
   for(int i=1;i<argc;i++)
        if(strcmp(argv[i],"-t")==0&&i+1<argc){
            nthreads=atoi(argv[i+1]);
            break;
        }
   
    solve(nthreads);
    
    return 0;
}
