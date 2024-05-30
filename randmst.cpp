#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <limits>
#include <ctime>
#include <functional>
#include <algorithm>
#include <chrono>
#include <random>
#include <cmath>

using namespace std;

static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<double> dis(0.0, 1.0);

// get a random weight between 0 and 1 for edges
double randWeight() {
    return dis(gen);
}

// priority queue class
class PriorityQueue {
private:
    vector<pair<int, double>> heap; // store pairs of vertex index and value
    vector<int> position; // track position of vertices in heap
    // binary heap
    void minHeapify(int i) {
        int smallest = i;
        int left = 2 * i + 1;
        int right = 2 * i + 2;
        // if left child is smaller than current smallest, update smallest
        if (left < heap.size() && heap[left].second < heap[smallest].second)
            smallest = left;
        // if right child smaller than current smallest, update smallest 
        if (right < heap.size() && heap[right].second < heap[smallest].second)
            smallest = right;
        // check if smallest element is the current, if not, then swap
        if (smallest != i) {
            // swap and update positions
            swap(position[heap[i].first], position[heap[smallest].first]);
            swap(heap[i], heap[smallest]);
            // recursively call minHeapify for the subtree
            minHeapify(smallest);
        }
    }

    void buildMinHeap() {
        for (int i = heap.size() / 2 - 1; i >= 0; --i) {
            minHeapify(i);
        }
    }

public:
    PriorityQueue(int n) : position(n, -1) {}
    // function to insert a new vertex and its key into the queue
    void insert(int vertex, double key) {
        position[vertex] = heap.size(); // update position of vertex
        heap.push_back({vertex, key}); // add element to end of heap
        int i = position[vertex]; // current index of inserted element

        // while not at root and parent's key is greater than the current key
        while (i != 0 && heap[(i - 1) / 2].second > heap[i].second) {
            swap(position[heap[i].first], position[heap[(i - 1) / 2].first]); // swap with parent, update position
            swap(heap[i], heap[(i - 1) / 2]);
            i = (i - 1) / 2;
        }
    }

    // return the verex with min key value 
    pair<int, double> deleteMin() {
        if (heap.empty()) return {-1, -1};
        // get min element and replace with last element in heap
        pair<int, double> root = heap[0];
        heap[0] = heap.back();
        position[heap[0].first] = 0; // update position of moved vertex
        heap.pop_back(); // remove last element
        position[root.first] = -1;

        if (!heap.empty()) minHeapify(0);

        return root; // return removed min element
    }

    // decrease key value
    void decreaseKey(int vertex, double key) {
        int i = position[vertex]; // i is the index of vertex in heap
        heap[i].second = key; // update key value
        // access second element or key, and check if current node's key is less than its parent's key 
        while (i != 0 && heap[(i - 1) / 2].second > heap[i].second) {
            swap(position[heap[i].first], position[heap[(i - 1) / 2].first]);
            // swap current node with its parent in heap
            swap(heap[i], heap[(i - 1) / 2]);
            //update current node's index, move up heap
            i = (i - 1) / 2;
        }
    }
    //check if vertex is currently in priority queue
    bool contains(int vertex) const {
        return position[vertex] != -1;
    }

    bool isEmpty() const {
        return heap.empty();
    }
};

// prim's algorithm to find the MST with n vertices
pair<double, bool> primMST(int n, function<double(int, int)> weightFunc, double threshold) {
    vector<double> key(n, numeric_limits<double>::max()); // store edge with min weight for each vertex
    vector<bool> inMST(n, false); // check if the vertex is in MST
    PriorityQueue pq(n); // object pq under PriorityQueue class to get vertex with the min key value

    // starting vertex
    int start = 0;
    key[start] = 0.0;
    pq.insert(start, 0.0); // insert the vertex in queue first 

    // loop
    while (!pq.isEmpty()) {
        int u = pq.deleteMin().first; // get the vertex with min key value
        inMST[u] = true; // mark vertex that it is in MST

        // update values and indices for adjacent vertices
        for (int v = 0; v < n; ++v) {
            if (u != v && !inMST[v]) {
                double weight = weightFunc(u, v); // get weight of edge from u to v
                // check if edge weight is below threshold
                if (weight <= threshold && weight < key[v]) { // if v is not in MST and weight of (u,v) is smaller than current key of v 
                    key[v] = weight; // update key value for v
                    if (pq.contains(v)) { // if v is already in queue, decrease key
                        pq.decreaseKey(v, weight);
                    } else {
                        pq.insert(v, weight); // otherwise insert v into priority queue
                    }
                }
            }
        }
    }
    // calculate total MST weight and check whether it's connected
    double totalWeight = 0.0;
    for (double w : key) {
        totalWeight += w;
    }

    // check if it includes all vertices
    bool isConnected = all_of(inMST.begin(), inMST.end(), [](bool v) { return v; });

    // primMST will output the total weight and connectivity
    return {totalWeight, isConnected};
}

// generate random points in 3D or 4D
vector<vector<double>> generateRandPoints(int n, int dimension, mt19937& gen, uniform_real_distribution<>& dis) {
    vector<vector<double>> points(n, vector<double>(dimension));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < dimension; ++j) {
            points[i][j] = dis(gen);
        }
    }
    return points;
}

// get Euclidean distance for 3D or 4D points
double euclideanDistance(const vector<double>& a, const vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += pow(a[i] - b[i], 2);
    }
    return sqrt(sum);
}

double handleThresholds(int n, function<double(int, int)> weightFunc, double initialThreshold, double biggerThreshold) {
    // find MST with the initial threshold
    auto initialResultPair = primMST(n, weightFunc, initialThreshold);
    // primMST outputs the total weight and whether it's connected
    double initialResult = initialResultPair.first; // total weight 
    bool isConnected = initialResultPair.second; // connectivity status

    // if graph isn't fully connected, try the bigger threshold
    if (!isConnected) {
        auto biggerResultPair = primMST(n, weightFunc, biggerThreshold);
        initialResult = biggerResultPair.first; // update initialResult with the new attempt
        isConnected = biggerResultPair.second; // update connection status
    }

    // return the result if connected, otherwise return infinity
    return isConnected ? initialResult : numeric_limits<double>::infinity();
}

int main(int argc, char* argv[]) {
    // checks if the 5 arguments are passed on command line
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << " 0 numpoints numtrials dimension" << endl;
        return 1;
    }

    // parse user's command line arguments
    int numpoints = stoi(argv[2]);
    int numtrials = stoi(argv[3]);
    int dimension = stoi(argv[4]);

    // initialize the threshold values based on the dimension
    double threshold = (7 - dimension) / (pow(numpoints, 0.5));
    double biggerThreshold = 80 * threshold;

    double averageWeight = 0.0;
    // track time
    auto startTime = chrono::high_resolution_clock::now();

    // loop runs numtrials times for a given dimension
    for (int i = 0; i < numtrials; ++i) {
        auto points = generateRandPoints(numpoints, dimension, gen, dis);
        // calculate edge weights
        auto weightFunc = [dimension, &points](int u, int v) {
            // if dimension is 0, generate random weight by calling randWeight()
            if (dimension == 0) {
                return randWeight();
            }
            return euclideanDistance(points[u], points[v]);
        };
        // calc MST weight. for each trial, add MST weight to avgWeight
        averageWeight += handleThresholds(numpoints, weightFunc, threshold, biggerThreshold);
    }

    // after the loop ends, divide the sum by number of trails 
    averageWeight /= numtrials;

    // stop timing and get time it took
    auto stopTime = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stopTime-startTime);
    
    // print output
    cout << averageWeight << " " << numpoints << " " << numtrials << " " << dimension << endl;
    // cout << "Time: " << duration.count() << " milliseconds " << endl;
    return 0;
}