
/// #######################################################################
/// A C++ program for Dijkstra's shortest path algorithm
/// The program is for adjacency matrix representation of the graph
/// You can choose between the N-th neighbors method or the epsilon method
/// Written by Francesco Tavanti, 2019 - 2020
/// #######################################################################


#include <stdio.h>
#include <iostream>
#include <array>
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <functional>
#include <list>
#include <utility>
#include <limits.h>
#include <fstream>
#include<cmath>




using namespace std;
using std::vector;

struct myclass {
    bool operator() (int i,int j) { return (i<j);}
} myobject;


//        ##################################
// An utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
//        ##################################

int minDistance(double dist[], bool sptSet[], int n)
{
    // Initialize min value
    int min = INT_MAX, min_index;
    
    for (int v = 0; v < n; v++){
        if (!sptSet[v]  && dist[v] <= min){
            min = dist[v];
            min_index = v;
        }
    }
    
    return min_index;
}


//        ##################################
// An utility function to print the constructed distance array
//        ##################################

void printSolution(double dist[], int n, int src){
    cout << "Vertex     Distance from Source"<<endl;
    double output_matrix [n][n];
    //    write output file
    ofstream ofs;
    ofs.open ("geodesic_out.dat", ofstream::out | ofstream::app);
    for (int i = 0; i < n; i++){
        cout << i+1 <<"  =====>   " << dist[i]<<endl;
        output_matrix[src][i]=dist[i];
        ofs << output_matrix[src][i] << ' ';
    }
    ofs << endl;
    ofs.close();
    
}


//        ##################################
//                   MAIN PROGRAM
//        ##################################

int main ()
{
    
    //    open file for input
    int m = 100;
//    m is the number of points in the matrix
    int iter=0,max_iter=m;
    double distances[m][m];
    ifstream myINfile;
    myINfile.open ("./geodesic.dat");
    cout<<endl;
    cout<<"Your input data"<<endl;
    for (int i = 0 ; i < m ; i++){
        for (int j = 0 ; j < m ; j++){
            myINfile >> distances[i][j];
        }
        
        cout<<endl;
    }
    myINfile.close();
    int n = m;
    
    ofstream ofs;
    ofs.open ("geodesic_out.dat", ofstream::out);
    
    int neigh_word;
    int nneigh;
    
    cout<<"****   INSTRUCTIONS   ****"<<endl<<endl;
    cout<<" This software reads the geodesic.dat file, which contains a matrix of distances "<<endl<<endl;
    cout<<" It will ask which method to use for the calculation of the path"<<endl;
    cout<<" In the Neighbors methods remember that not all nodes are connected each other so"<<endl;
    cout<<" do not chose a too high number of neighbors (usually less than max number of nodes/2)"<<endl<<endl;
    cout<<" In the Epsilon method if you chose a too small value the path will finish when there will be no available nodes"<<endl<<endl;
    cout<<" In output the file geodesic_out.dat contains the path obtained and can be plotted"<<endl<<"****    ****"<<endl;
    
    cout <<"If you want the Neighbors method insert 1 "<<endl;
    cout <<"If you want the Epsilon method insert 0 "<<endl;
    cin >>neigh_word;
    
    
    double distance_neighbors;
    
    //    define all distances to infinite
    double h[n][n];
    bool indices[n][n];
   
    int new_parent[n];
    int old_parent[n];
    
    //  set all values in h to INF
    double a = numeric_limits<double>::infinity();
    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            h[i][j]=1*a;
        }
    }
    
    if(neigh_word==1){
        //        ##################################
        //        This performs the Neighbors method
        //        ##################################
        cout<<endl;
        cout<<"Insert the number of nearest neighbors: "<<endl;
        cin >> nneigh;
        
        if(nneigh>(m/2))
        {
            cout<<"The number of nearest neighbors chosen "<<nneigh<<" is greater than the number of nodes "<<m<<endl;
            exit(EXIT_FAILURE);
        }
        cout <<"Performing the Neighbors method"<<endl;
        //    variables initialization
        int k=0;
        cout<<endl;
        
        
        // vectors to sort and to store the neighbors
        vector<vector<pair<double,int> > > neighbors;
        neighbors.resize(n);
        for (int i=0; i<n; i++) {
            for (int j=0; j<m; j++) {
                neighbors[i].push_back(make_pair(distances[i][j],k));
                k++;
            }
            sort(neighbors[i].begin(),neighbors[i].begin()+m);
        }
        
        for (int i=0;i<m;i++){
            for (int j=0;j<m;j++){
                distance_neighbors=neighbors[i][j].first;
                k=(neighbors[i][j].second)-n*i;
      
                //    -n*i-1 rescale the neighbors for the size of the data, so each row starts from 0
                if(j<=nneigh){
                    h[i][k]=distances[i][k];
                    indices[i][k]=false;
                    
                }
                else{
                    h[i][k]=a;
                    indices[i][k]=true;
                }
               
            }
         
        }
        
    }
    //        ##################################
    //        End of the neighbors method
    //        ##################################
    
    else{
        
        //        ##################################
        //        This performs the epsilon method
        //        ##################################
        double epsilon_value;
        cout<<endl;
        cout <<"Insert a value for epsilon: "<<endl;
        cin >>epsilon_value;
        if(epsilon_value<=0)
        {
            cout<<"Your MAXCONNECT value is equal or lower than 0 and this has no physical meaning. Change the value to a positive value"<<endl;;
        exit(EXIT_FAILURE);
        }
        cout <<"Performing the Epsilon method"<<endl;
        cout<<endl;
        cout<< "Matrix of distances "<<endl;
        
        double dist;
        for (int i=0;i<m;i++){
            
            for (int j=0;j<m;j++){
                dist=distances[i][j];
              
                if(dist<=epsilon_value){
                    
                    h[i][j]=distances[i][j];
                    indices[i][j]=false;
                }
                else{
                    h[i][j]=a;
                    indices[i][j]=true;
                }
                cout<< h[i][j]<<"  ";
            }
            cout << endl;
        }
    }
    //        ##################################
    //        End of the epsilon method
    //        ##################################
    
    
    //        ##################################
    // Function that implements Dijkstra's single source shortest path algorithm
    // for a graph represented using adjacency matrix representation
    //        ##################################
    
    int j=0;
    int src=j;
    double dist[n];
    // The output array.  dist[i] will hold the shortest
    // distance from j to i
    
    bool sptSet[n];
    // sptSet[i] will be true if vertex i is included in shortest
    // path tree or shortest distance from j to i is finalized
    // Initialize all distances as INFINITE and stpSet[] as false
    for (int i = 0 ; i < n ; i++){
        dist[i] = a;
        sptSet[i] = false;
        
    }
    
    // Distance of source vertex from itself is always 0
    dist[src] = 0;
    // Find shortest path for all vertices
    for (int count = 0; count < n; count++)
    {
        // Pick the minimum distance vertex from the set of vertices not
        // yet processed. u is always equal to j in the first iteration.
        int u = minDistance(dist, sptSet, n);
        // Mark the picked vertex as processed
        sptSet[u] = true;
        
        // Update dist value of the adjacent vertices of the picked vertex.
        for (int v = 0; v < n; v++){
            
            // Update dist[v] only if is not in sptSet, there is an edge from
            // u to v, and total weight of path from src to  v through u is
            // smaller than current value of dist[v]
            if (!sptSet[v] && h[u][v] && dist[u] != INT_MAX
                && dist[u]+h[u][v] < dist[v]){
                
                dist[v] = dist[u] + h[u][v];
                old_parent[v]=u;
                new_parent[v]=v;

                if(old_parent[v]==v)
                {
                    break;
                }
                else
                {
                    iter++;
                }
                
                cout << "vertex "<< old_parent[v] << "   to  vertex " << new_parent[v] << "     distance --> " << dist[v]<<" ITERATION "<<iter<< endl;
                ofs<<old_parent[v]<<" "<<new_parent[v]<<" "<<iter<<endl;
                break;
                
            }
            
        }
    
    }

    ofs << endl;
    ofs.close();
    return 0;
    
}

