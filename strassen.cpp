#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <random>

using namespace std;

typedef vector<vector<int>> Matrix; 

// function to pad matrix with zeros
void padMatrix(Matrix& matrix) {
    int originalSize = matrix.size();
    int newSize = originalSize + 1;
    for (int i = 0; i < originalSize; ++i) {
        matrix[i].push_back(0);
    }
    matrix.push_back(vector<int>(newSize, 0));
}

// function used to read matrix from the input file
Matrix readMatrix(std::ifstream& inputFile, int dimension) {
    // create empty matrix
    Matrix matrix(dimension, std::vector<int>(dimension));
    // fill in matrix with values read from the input file
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            if (!(inputFile >> matrix[i][j])) {
                std::cout << "Error reading matrix from file. Check the file format and dimensions." << std::endl;
                return {};  // Return an empty matrix in case of error
            }
        }
    }
    return matrix;
}

// output diagonal elements of a matrix
void outputDiagonal(const Matrix& matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        std::cout << matrix[i][i] << std::endl;  // output the diagonal element at position (i, i)
    }
}
// multiply 2 matrices using conventional method
Matrix conventionalMultiply(const Matrix& A, const Matrix& B) {
    size_t n = A.size(); // get size of matrices
    Matrix C(n, vector<int>(n, 0)); // empty result matrix
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < n; k++) {
            for (size_t j = 0; j < n; j++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}
// add two matrices
Matrix add(const Matrix& A, const Matrix& B) {
    size_t n = A.size();
    Matrix result(n, vector<int>(n, 0));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}
// subtract corresponding elements of second from first matrix
Matrix subtract(const Matrix& A, const Matrix& B) {
    size_t n = A.size(); 
    Matrix result(n, vector<int>(n, 0));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}

// combine results of Strassen's sub problems
Matrix combine(const Matrix& C11, const Matrix& C12, const Matrix& C21, const Matrix& C22) {
    int n = C11.size();
    Matrix C(2 * n, vector<int>(2 * n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = C11[i][j];
            C[i][j + n] = C12[i][j];
            C[i + n][j] = C21[i][j];
            C[i + n][j + n] = C22[i][j];
        }
    }
    return C;
}
// splits the matrix into four submatrices
void split(const Matrix& original, Matrix& a11, Matrix& a12, Matrix& a21, Matrix& a22) {
    int newSize = original.size() / 2;
    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            a11[i][j] = original[i][j];
            a12[i][j] = original[i][j + newSize];
            a21[i][j] = original[i + newSize][j];
            a22[i][j] = original[i + newSize][j + newSize];
        }
    }
}

void removePadding(Matrix& matrix, int originalSize) {
    // resize each row to the original size
    for (auto& row : matrix) {
        row.resize(originalSize);
    }
    // resize the matrix itself to have the original number of rows
    matrix.resize(originalSize);
}

// Multiply two matrices using Strassen's algorithm, adjusting for a dynamic crossover point.
Matrix strassenMultiply(Matrix A, Matrix B, int originalSize) {
    int currentSize = A.size();

    // If the matrix is 1x1, multiply the single elements directly.
    if (currentSize == 1) {
        return Matrix(1, vector<int>(1, A[0][0] * B[0][0]));
    }

    // Determine the crossover point dynamically based on the original size of the matrix.
    // For matrices with an even original size, the algorithm will switch to conventional multiplication
    // after one iteration of Strassen's algorithm. For matrices with an odd original size, the switch
    // occurs when the size is reduced to (originalSize + 1)/2
    int crossoverPoint = (originalSize % 2 == 0) ? (originalSize / 2) : ((originalSize + 1) / 2);

    // Use conventional multiplication if the current size of the matrix is at or below the crossover point.
    if (currentSize <= crossoverPoint) {
        return conventionalMultiply(A, B);
    }

    // Split the current matrices into four submatrices each, preparing for recursive multiplication.
    int newSize = currentSize / 2;
    Matrix a11(newSize, vector<int>(newSize)), a12(newSize, vector<int>(newSize)),
           a21(newSize, vector<int>(newSize)), a22(newSize, vector<int>(newSize)),
           b11(newSize, vector<int>(newSize)), b12(newSize, vector<int>(newSize)),
           b21(newSize, vector<int>(newSize)), b22(newSize, vector<int>(newSize));

    // Splitting the matrices A and B into submatrices a11, a12, a21, a22 and b11, b12, b21, b22 respectively.
    split(A, a11, a12, a21, a22);
    split(B, b11, b12, b21, b22);

    // Recursively multiply submatrices using Strassen's algorithm, calculating seven products.
    Matrix P1 = strassenMultiply(add(a11, a22), add(b11, b22), originalSize);
    Matrix P2 = strassenMultiply(add(a21, a22), b11, originalSize);
    Matrix P3 = strassenMultiply(a11, subtract(b12, b22), originalSize);
    Matrix P4 = strassenMultiply(a22, subtract(b21, b11), originalSize);
    Matrix P5 = strassenMultiply(add(a11, a12), b22, originalSize);
    Matrix P6 = strassenMultiply(subtract(a21, a11), add(b11, b12), originalSize);
    Matrix P7 = strassenMultiply(subtract(a12, a22), add(b21, b22), originalSize);

    // Combine the products to form the resulting matrix, according to Strassen's formula.
    Matrix c11 = add(subtract(add(P1, P4), P5), P7);
    Matrix c12 = add(P3, P5);
    Matrix c21 = add(P2, P4);
    Matrix c22 = add(subtract(add(P1, P3), P2), P6);

    // Combine the four submatrices into one final matrix and return it.
    return combine(c11, c12, c21, c22);
}


Matrix strassen(Matrix A, Matrix B) {
    int originalSize = A.size(); // This is the original size of the matrices

    // Handle odd-sized matrices by padding them
    bool odd = originalSize % 2 != 0;
    if (odd) {
        padMatrix(A);
        padMatrix(B);
    }

    // Call strassenMultiply with the original size of the matrices
    Matrix result = strassenMultiply(A, B, originalSize);

    // If matrices were padded, remove the padding from the result
    if (odd) {
        removePadding(result, originalSize);
    }

    return result;
}


std::chrono::duration<double, std::milli> trackTime(const Matrix& a, const Matrix& b, bool useStrassen) {
    // start the timer
    auto start = std::chrono::high_resolution_clock::now();
    // useStrassen is a boolean, if true we use strassen and false then we use conventional way
    if (useStrassen){
        strassen(a,b);
    } else{
        conventionalMultiply(a,b);
    }
    // end the timer
    auto end = std::chrono::high_resolution_clock::now();
    return end - start;
}

Matrix generateRandomMatrix(int n) {
    Matrix matrix(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = rand() % 2; // random values between 0 and 1
        }
    }
    return matrix;
}

// find the smallest odd/even size of a matrix in which Strassen's algorithm is faster than conventioanl
// if strassen's algorithm is not faster, we increment i by 2 (to the next odd/even size) until we find the crossover point
int findCrossover(bool startWithOdd){
    // do true to start at odd-sized matrix; false to start at even-sized matrix
    int start = startWithOdd ? 1 : 2;
    
    for (int i = start; ; i += 2) { // increment by 2 to either stay odd or even
        Matrix A = generateRandomMatrix(i);
        Matrix B = generateRandomMatrix(i);
        auto time1 = trackTime(A, B, false); // conventional multiplication
        printf("Time for conventional: %f, n: %d\n", time1.count(), i);
        auto time2 = trackTime(A, B, true);  // Strassen's algorithm
        printf("Time for Strassen's: %f, n: %d\n", time2.count(), i);
        
        // if Strassen's algorithm is faster, return the current size as the crossover point
        if (time1.count() > time2.count()) {
            return i;
        }
    }
}


// random adjacency matrix for a graph
vector<vector<int>> generateRandomGraph(int n, double p) {
    vector<vector<int>> graph(n, vector<int>(n, 0)); // initialize n by n matrix to 0
    // random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution dist(p); // boolean distribution with probability p

    // loop through the upper triangle of the matrix 
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (dist(gen)) { // if generated random boolean is true w prob p
                graph[i][j] = 1; // set edge to 1 
                graph[j][i] = 1; // set edge to 1 for undirected graph
            }
        }
    }
    return graph;
}
// count triangles using strassen's
long long countTriangles(const vector<vector<int>>& graph) {
    auto A2 = strassen(graph, graph); // square adjacency matrix
    auto A3 = strassen(A2, graph); // cube it

    long long triangles = 0;
    // sum the diagonals of A3. These are each triangles touching vertex i
    for (int i = 0; i < A3.size(); ++i) {
        triangles += A3[i][i];
    }
    return triangles / 6; // divide by 6 since each triangle is counted 6 times
}

double nChooseK(int n, int k) {
    double result = 1;
    for (int i = 1; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

// main function
int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <flag> <dimension> <inputfile>" << std::endl;
        return 1;
    }
    // command line inputs
    int flag = std::stoi(argv[1]);
    int dimension = std::stoi(argv[2]);
    std::string inputFilePath = argv[3];
    std::ifstream inputFile(inputFilePath);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file" << std::endl;
        return 1;
    }

    // print crossover value!!
    // int crossoverPoint1 = findCrossover(false);
    // std::cout << "Crossover point: " << crossoverPoint1 << std::endl;

    // read matrices from the input file
    Matrix A = readMatrix(inputFile, dimension);
    Matrix B = readMatrix(inputFile, dimension);
    inputFile.close();

    // perform matrix multiplication using conventional way
    Matrix result = conventionalMultiply(A, B);

    // get diagonal elements of the result matrix
    outputDiagonal(result);

    // triangle stuff
    /* int n = 1024;
    vector<double> ps = {0.01, 0.02, 0.03, 0.04, 0.05};

    cout << "Probability | Triangles (1 Trial) | Triangles (5 Trials) | Expected Triangles\n";
    cout << "-----------------------------------------------------------------------------\n";

    for (double p : ps) {
        auto graph = generateRandomGraph(n, p);
        long long triangles1 = countTriangles(graph); // 1 trial

        long long sumTriangles = 0;
        for (int trial = 0; trial < 5; ++trial) { // do 5 trials
            auto graph = generateRandomGraph(n, p);
            sumTriangles += countTriangles(graph); 
        }
        long long avgTriangles = sumTriangles / 5; // average the count over the 5 trials

        long long expectedTriangles = nChooseK(n, 3) * p * p * p;
        cout << "p = " << p << " | " << triangles1 << " | " << avgTriangles << " | " << expectedTriangles << endl;
    }*/
    return 0;
}