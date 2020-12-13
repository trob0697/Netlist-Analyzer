/*
 *  By: Tahreak Robinson
 *  Date: 26 November 2020
 *  Reference 1: https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.571.2560&rep=rep1&type=pdf
 *  Reference 2: https://en.wikipedia.org/wiki/Row_echelon_form
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
using namespace std;

vector< vector<string> > components;
int numNodes = 0, numBranches = 0, size = 0;
vector< vector<double> > incidenceMatrix;
vector< vector<double> > voltageCoefMatrix;
vector< vector<double> > currentCoefMatrix;
vector< vector<double> > sparseMatrix;
vector<double> circuitSolution;
bool isSingularMatrix = false;

template <typename T>
void debugPrintMatrix(vector< vector<T> > v){
    for(int i = 0; i < v.size(); i++){
        for(int j = 0; j < v[0].size(); j++){
            cout << v[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

template <typename T>
void debugPrintArray(vector<T> v){
    for(int i = 0; i < v.size(); i++){
        cout << v[i] << "\t";
    }
    cout << endl;
}

void clearNetList(){
    components.clear();
    numNodes = 0, numBranches = 0, size = 0;
    incidenceMatrix.clear();
    voltageCoefMatrix.clear();
    currentCoefMatrix.clear();
    sparseMatrix.clear();
    circuitSolution.clear();
    isSingularMatrix = false;
}

int readNewNetlist(){
    cout << "Enter the file name: ";
    string fileName;
    cin >> fileName;
    fileName = "./netlists/" + fileName;

    ifstream myFile(fileName);
    if(myFile.is_open()){
        string line;
        while(getline(myFile, line)){
            vector<string> comp;
            string w = "";
            for(int i = 0; i <= line.length(); i++){
                if(line[i] == ' ' || i == line.length()){
                    comp.push_back(w);
                    w = "";
                }
                else{
                    w += line[i];
                }
            }
            components.push_back(comp);
        }

        myFile.close();
        return 1;
    }
    else{
        cout << "Unable to open file" << endl;
        return 0;
    }
}

void countNumBranchesAndNodes(){
    numBranches = components.size();
    int sourceNodeNum, destinationNodeNum, largestNodeNum;
    for(int i = 0; i < components.size(); i++){
        sourceNodeNum = stoi(components[i][1]);
        destinationNodeNum = stoi(components[i][2]);
        largestNodeNum = max(sourceNodeNum, destinationNodeNum);
        numNodes = max(numNodes, largestNodeNum);
    }
    size = numNodes + numBranches*2;
}

void createIncidenceMatrix(){
    vector< vector<double> > vec(numNodes, vector<double> (numBranches, 0));
    incidenceMatrix = vec;

    for(int branch = 0; branch < numBranches; branch++){
        for(int j = 1; j <= 2; j++){
            int node = stoi(components[branch][j]);
            int val;
            if(node != 0){
                if(j == 1){
                    val = 1;
                }
                else{
                    val = -1;
                }
                incidenceMatrix[node-1][branch] = val;
            }
        }
    }
}

void createVoltageCoefMartix(){
    vector< vector<double> > vec(numBranches, vector<double> (numBranches, 0));
    voltageCoefMatrix = vec;
    for(int i = 0; i < numBranches; i++){
        voltageCoefMatrix[i][i] = 1;
    }
}

void createCurrentCoefMatrix(){
    vector< vector<double> > vec(numBranches, vector<double> (numBranches, 0));
    currentCoefMatrix = vec;
    for(int i = 0; i < numBranches; i++){
        char componentType = components[i][0].at(0);
        if(componentType == 'V'){
            currentCoefMatrix[i][i] = 0;
        }
        else if(componentType == 'R'){
            currentCoefMatrix[i][i] = stod(components[i][3]) * (-1);
        }
    }
}

void createSparseMatrix(){
    vector< vector<double> > vec(size, vector<double> (size, 0));
    sparseMatrix = vec;
    
    int i, j, k, l;
    //Input Incidence Matrix
    k = 0;
    for(i = 0; i < numNodes; i++){
        l = 0;
        for(j = numNodes+numBranches; j < size; j++){
            sparseMatrix[i][j] = incidenceMatrix[k][l];
            l++;
        }
        k++;
    }
    
    //Input Negative Transposed Incidence Matrix
    k = 0;
    for(i = numNodes; i < numNodes+numBranches; i++){
        l = 0;
        for(j = 0; j < numNodes; j++){
            sparseMatrix[i][j] = incidenceMatrix[l][k] * (-1);
            l++;
        }
        k++;
    }

    //Input Identity Matrix
    for(i = numNodes; i < numNodes+numBranches; i++){
        sparseMatrix[i][i] = 1;
    }

    //Input Voltage Coefficient Matrix
    k = 0;
    for(i = numNodes+numBranches; i < size; i++){
        l = 0;
        for(j = numNodes; j < numNodes+numBranches; j++){
            sparseMatrix[i][j] = voltageCoefMatrix[k][l];
            l++;
        }
        k++;
    }

    //Input Current Coefficient Matrix
    k = 0;
    for(i = numNodes+numBranches; i < size; i++){
        l = 0;
        for(j = numNodes+numBranches; j < size; j++){
            sparseMatrix[i][j] = currentCoefMatrix[k][l];
            l++;
        }
        k++;
    }
}

void appendInputToSpareMatrix(){
    int i = 0, j =0;
    while(i < sparseMatrix.size()){
        int val = 0;
        if(i >= numNodes+numBranches){
            if(components[j][0].at(0) == 'V'){
                val = stod(components[j][3]);
            }
            j++;
        }
        sparseMatrix[i].push_back(val);
        i++;
    }
}

/******************************************************************/
/********************** Guassian Elimination **********************/
/******************************************************************/

// Function to convert matrix to R.R.E.F.
int rref(){
    for(int i = 0; i < size; i++){
        // Find the row with the max value at the i-th column
        int i_max = i;
        int i_max_val = sparseMatrix[i][i];
        for(int k = i; k < size; k++){
            int abs_cur_val = abs(sparseMatrix[k][i]);
            if(abs_cur_val > i_max_val){
                i_max = k;
                i_max_val = abs_cur_val;
            }
        }

        //If the max value is zero, the matrix is singular
        if(i_max_val == 0){
            return i;
        }

        // If the row with the max value isnt on top, swap with the top row
        if(i_max != i){
            vector<double> temp = sparseMatrix[i];
            sparseMatrix[i] = sparseMatrix[i_max];
            sparseMatrix[i_max] = temp;
        }

        // Reduce max row
        double d = sparseMatrix[i][i];
        for(int I = i; I <= size; I++){
            sparseMatrix[i][I] /= d;
        }

        // Reduce other rows
        for(int I = i+1; I < size; I++){
            if(sparseMatrix[I][i] == 0){
                continue;
            }
            double temp = sparseMatrix[I][i];
            for(int J = i; J <= size; J++){
                sparseMatrix[I][J] += (-1)*temp*sparseMatrix[i][J];
            }
        }
    } 
    return -1; 
}

//Function to calculate the values of the unknowns
void backSubstitutioh(){ 
    // Initialize solution array
    vector<double> vec(size, 10);
    circuitSolution = vec;

    // Start calculating from last equation up to the first
    for(int i = size-1; i >= 0; i--){ 
        // Start with the RHS of the equation
        circuitSolution[i] = sparseMatrix[i][size]; 
  
        // Initialize j to i+1 since matrix is upper triangular
        for(int j = i+1; j < size; j++){ 
            /* subtract all the lhs values 
             * except the coefficient of the variable 
             * whose value is being calculated */
            circuitSolution[i] -= sparseMatrix[i][j]*circuitSolution[j]; 
        } 

        // Divide the RHS by the coefficient of the unknown being calculated
        circuitSolution[i] = circuitSolution[i]/sparseMatrix[i][i]; 
    } 
}

void guassianElimination(){
    // Reduction into R.E.F.
    int singular_flag = rref();
  
    // If matrix is singular
    if(singular_flag != -1){ 
        isSingularMatrix = true;
        // If the RHS of equation corresponding to zero row  is 0, * system has infinitely many solutions, else inconsistent
        if(sparseMatrix[singular_flag][size]) 
            cout << "Singular Matrix: Inconsistent System" << endl; 
        else
            cout << "Singular Matrix: May have infinitely many solutions." << endl; 
  
        return; 
    } 
    else{
        // Get solution to system
        isSingularMatrix = false;
        backSubstitutioh();
    }
}

/******************************************************************/
/******************************************************************/

void printCurrentValues(){
    for(int i = numNodes+numBranches; i < size; i++){
        cout << "I" << i-(numNodes+numBranches)+1 << " = " << circuitSolution[i] << "A" << endl;
    }
}

void printVoltageValues(){
    cout << "E0 = 0V" << endl;
    for(int i = 0; i < numNodes; i++){
        cout << "E" << i+1 << " = " << circuitSolution[i] << "V" << endl;
    }
    for(int i = numNodes; i < numNodes+numBranches; i++){
        cout << "V" << i-numNodes+1 << " = " << circuitSolution[i] << "V" << endl;
    }
}

int main(){
    bool running = true;
    while(running){
        cout << "Select one of the following options" << endl;
        cout << "A. Read  new netlist" << endl;
        cout << "B. Compute current values for current netlist" << endl;
        cout << "C. Compute voltage values for current netlist" << endl;
        cout << "D. Exit" << endl;
        cout << "> ";

        char option;
        cin >> option;
        cout << endl;

        switch(option){
            case 'A':
                cout << "You have selected \"Read new netlist\"" << endl;
                clearNetList();
                if(readNewNetlist()){
                    countNumBranchesAndNodes();
                    createIncidenceMatrix();
                    createVoltageCoefMartix();
                    createCurrentCoefMatrix();
                    createSparseMatrix();
                    appendInputToSpareMatrix();
                    guassianElimination();
                }
                break;
            case 'B':
                cout << "You have selected \"Compute current values for current netlist\"" << endl;
                if(isSingularMatrix){
                    cout << "Please select a valid netlist" << endl;
                }
                else{
                    printCurrentValues();
                }
                break;
            case 'C':
                cout << "You have selected \"Compute voltage values for current netlist\"" << endl;
                if(isSingularMatrix){
                    cout << "Please read a valid netlist" << endl;
                }
                else{
                    printVoltageValues();
                }
                break;
            case 'D':
                cout << "Exiting program" << endl;{}
                running = false;
                break;
            default:
                cout << "Invalid input" << endl;
                break;
        }
        cout << endl;
    }

    return 0;
}