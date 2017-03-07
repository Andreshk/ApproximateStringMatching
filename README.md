# Levenshtein distance (edit distance)
The [`levenshtein`](https://github.com/Andreshk/ApproximateStringMatching/blob/master/Source.cpp#L17) function calculates the well-known and established Levenshtein edit distance between two strings. It is a simple dynamic programming algorithm that uses a 2D matrix for memoization and, for strings of length `m` and `n`, has `O(mn)` time complexity and `O(min(m,n))` space complexity. The version presented here outputs the full matrix and recovers the entire edit process, visualizing it for user convenience (therefore using `O(mn)` space).

Transposition of adjacent characters can be allowed with a boolean parameter, transforming the calculated distance into Optimal String Alignment distance. It should be noted that this differs from [Damerau-Levenshtein](https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance) distance and is not a metric, whereas the Levenshtein distance is. Nevertheless, the algorithm calculating Levenshtein distance can be easily modified to allow transpositions. 

# Approximate String Matching

Almost all approximate string matching algorithms are based on the Levenshtein distance algorithm, with a small modification to the DP matrix. There are 4 algorithms presented here, all of which only return the # of "matches" found with a limit of `k` errors. It is well-established to use `m` for the pattern length (the string to search) and `n` for the text length (the string in which we search). It is assumed that `m` is, often by magnitudes, smaller than `n`.

## Classic dynamic programming

The full dynamic programming algorithm calculates the entire table column by column, only keeping the current and the previous column in memory. It thus takes `О(mn)` time and uses `O(m)` space.

## Dynamic programming with cut-off heuristic

This algorithm does not calculate every single column, if it judges that after some point the values in this column will not affect the result. It is postulated and it can be shown that the average running time is `O(kn)`, although in the worst case it remains `O(mn)`. The space complexity is still `O(m)`.

## Diagonal transitions algorithm

An important property of the matrix are its non-decreasing diagonals (parallel to the main diagonal). This algorithm manages to compute only the places along each diagonal at which the values increase (for example, 0 0 0 **1** 1...). The time complexity is again `O(kn)` on average, `O(mn)` in the worst case, but the space complexity is `O(k)`.

## Bitwise parallelized matrix

This unusual algorithms exploits some finer properties of the matrix and manages the store the state of each column in only a couple of binary vectors. The calculations for the next column can be simulated with bitwise operations over these vectors, making them suitable for storage in machine words. For simplicity they are limited to 64 bits, meaning that `m` can not exceed 64. If we take the size of the alphabet `b` into consideration, we have `O(bm + n)` time and `O(bm)` complexity. In practise `b` is most often limited by a constant, and with `m` up to 64 we can think of the complexities as `O(n)` and `O(1)` respectively.
