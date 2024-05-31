# prog1
Overview of My MST Algorithm:
To find the minimum spanning tree (MST) for graphs in varying dimensions, I implemented Prim's algorithm, chosen for its efficiency in graphs with a high edge-to-vertex ratio. I adapted it by incorporating a binary heap for the priority queue, drawing on pseudocode from lecture materials. My primary functions within this algorithm include insert, for adding new vertex-edge pairs, decreaseKey for updating edge weights, and deleteMin for extracting the minimum edge.

Optimizing the Process:
Initially, I worked with the basic structure of the algorithm without optimizations. To enhance performance, I introduced a thresholding technique to prune unnecessary edges, thereby reducing computational load. The threshold was based on the number of vertices (numpoints) and the graph's dimensionality, set as $\frac{7 - \text{dimension}}{\sqrt{numpoints}}$, with adjustments through empirical testing. If the first run of the algorithm does not yield a connected MST, it retries with an increased threshold.

This approach not only streamlined the computation by focusing on potentially useful edges but also allowed dynamic adaptation to different graph sizes and dimensions. My further improvement plan includes exploring a randomized approach to further reduce the computation area by focusing on a subset of possible connections.

To run: After git cloning the repo with the makefile and python file, and cd'ing, execute python3 prog1.py in the command line.

# prog2
In my exploration to determine the cross-over point between conventional matrix multiplication and Strassen's algorithm, I conducted experiments to identify at which matrix size Strassen's becomes more efficient. My method, implemented in a function called trackTime, dynamically adjusts for odd and even-sized matrices and tracks performance, pinpointing when Strassen's algorithm starts outperforming the conventional method.

For odd-sized matrices, the most consistent cross-over point identified was at size 129, while for even-sized matrices, it was at size 64. These points varied slightly with each trial due to factors like matrix randomness and background computing tasks, yet they remained around these values consistently.

This discovery was higher than my initial theoretical calculations, which didn't account for the overhead of memory allocation and the copying of matrix values—operations that Strassen's algorithm performs frequently. In practice, these operations are significant, especially for smaller matrices where the overhead can outweigh the algorithm's theoretical efficiency.

Additionally, I noted that odd-sized matrices required more padding (to make dimensions even for Strassen's algorithm), which further increases the computation time due to extra recursive calls and memory allocation.

I also optimized the conventional matrix multiplication method by reordering the looping variables to improve cache memory utilization, enhancing performance significantly.

By dynamically choosing between Strassen's and conventional methods based on matrix size and observing these cross-over points, I aimed to optimize overall algorithm performance. If I had more resources, I'd further optimize by pre-allocating memory for matrix results, reducing the overhead from recursive memory allocation in Strassen's algorithm.

This exploration underlines that while theoretical models provide a baseline, practical implementation often reveals additional complexities that must be managed to truly optimize performance.

# prog3

In my analysis of seven main algorithms to find the optimal partitioning residue, I considered their effectiveness and efficiency based on the number of iterations, $M$, where each algorithm behaves differently. Here’s a summary from my findings and comparisons:

KK Algorithm: As a deterministic approach, KK consistently produces near-optimal solutions quickly due to its $\mathcal{O}(n\log n)$ complexity. However, it often fails to achieve the true minimum because it forces the largest two elements into opposing subsets. From my data, KK generally outputs residues around $2 \times 10^5$.

Repeated Random: This randomized method shows varying outcomes for the same input. It’s relatively slow, requiring many iterations to approach a good residue, typically achieving residues around $2 \times 10^9$ for 25,000 iterations. The method would need an impractically high number of iterations, estimated at $2.5 \times 10^8$, to match the efficiency of the KK algorithm.

Hill Climbing: This heuristic method updates residues efficiently in $\mathcal{O}(1)$ time per iteration but tends to get stuck at local minima, making it less reliable for higher iterations.

Simulated Annealing: Similar to hill climbing but allows transitions to worse states, potentially escaping bad local minima. Although it requires more iterations to match hill climbing’s results, it outperforms in avoiding poor local solutions.

Pre-partitioned Repeated Random: Incorporating random pre-partitioning with the efficient KK algorithm in each iteration slightly lengthens the runtime compared to standard repeated random but shows potential for superior results due to the combination of randomization and deterministic reduction.

Pre-partitioned Hill Climbing and Simulated Annealing: These methods apply their respective strategies to pre-partitioned sets, sharing similar time complexities with their non-pre-partitioned counterparts but showing varied effectiveness based on the initial partitions.

Overall, the KK algorithm provides a reliable and quick approximation to the optimal solution but lacks the potential to consistently find the true minimum. In contrast, randomized approaches, especially those involving pre-partitioning, offer a higher chance of approaching or surpassing the KK algorithm’s performance given sufficient iterations. Pre-partitioned methods particularly stand out by effectively combining randomization with deterministic strategies to achieve significantly better outcomes.

Using the KK output as a starting point for other algorithms can enhance their efficiency, especially in scenarios where a good initial approximation leads to quicker convergence on optimal or near-optimal solutions. This strategy is particularly beneficial for hill climbing, which can benefit from starting near a local minimum but might otherwise get stuck without a strong starting position. In summary, while each method has its merits and limitations, the choice of algorithm and its configuration can greatly influence the effectiveness and efficiency of finding the best partitioning solution.
