#include<iostream>
#include<iomanip>
#include<vector>
#include<string> // std::to_string
#include<string_view>
#include<type_traits>

template<typename T>
T min(T a, T b) { return a < b ? a : b; }
template<typename T>
T min(T a, T b, T c) { return a < b ? (a < c ? a : c) : (b < c ? b : c); }
template<typename T>
T max(T a, T b, T c) { return a > b ? (a > c ? a : c) : (b > c ? b : c); }

// Calculates classic Levenshtein distance between two strings. By allowing
// neighbouring character transpositions it turns into Optimal String Alignment distance.
// (Note: this differs from Damerau-Levenshtein distance, and unlike Levenshtein is not a metric.)
void levenshteinDemo(std::string_view p, std::string_view t, bool withTranspositions = false) {
    const size_t m = p.length() + 1, n = t.length() + 1;
    // simulating a matrix with table[i][j] <=> array[i*n + j]
    std::vector<size_t> table(m*n);
#define _table(i,j) table[(i)*n + (j)]
    for (size_t i = 0; i < m; i++) _table(i, 0) = i;
    for (size_t j = 0; j < n; j++) _table(0, j) = j;

    for (size_t i = 1; i < m; i++)
        for (size_t j = 1; j < n; j++) {
            _table(i, j) = min(_table(i - 1, j) + 1,
                               _table(i, j - 1) + 1,
                               _table(i - 1, j - 1) + 1 - (p[i - 1] == t[j - 1]));
            // View transposition as an option only if the last 2 characters are the same
            if (withTranspositions && i > 1 && j > 1 && p[i - 1] == t[j - 2] && p[i - 2] == t[j - 1])
                _table(i, j) = min(_table(i, j), _table(i - 2, j - 2) + 1);
        }

    // return _table(m - 1, n - 1);
    // Pretty printing the entire matrix
    const size_t numw = std::to_string(m > n ? m : n).length();
    std::cout << std::string(numw + 2, ' ');
    for (size_t j = 0; j < n - 1; j++)
        std::cout << " " << std::setw(numw) << t[j];
    std::cout << "\n";
    for (size_t i = 0; i < m; i++) {
        std::cout << ((i == 0) ? ' ' : p[i - 1]);
        for (size_t j = 0; j < n; j++)
            std::cout << " " << std::setw(numw) << _table(i, j);
        std::cout << "\n";
    }
    std::cout << "Distance = " << _table(m - 1, n - 1) << "\n";

    // Reverting the entire edit process, a single operation at a time
    std::string p1, t1;
    p1.reserve(m + n);
    t1.reserve(m + n);
    // Starting from the bottom right corner, at each step we recover
    // the last performed operation by checking the current value
    // from which of its neighbouring values has been obtained.
    size_t i = m - 1, j = n - 1;
    while (i && j) {
        if (withTranspositions && i > 1 && j > 1
            && _table(i, j) == _table(i - 2, j - 2) + 1
            && p[i - 1] == t[i - 2] && p[i - 2] == t[j - 1])
        { // The last operation has been a successful (and allowed) transposition
            p1.push_back('}'); p1.push_back(p[i - 1]); p1.push_back(p[i - 2]); p1.push_back('{');
            t1.push_back('}'); t1.push_back(t[j - 1]); t1.push_back(t[j - 2]); t1.push_back('{');
            i -= 2; j -= 2;
        } else if (_table(i, j) == _table(i - 1, j - 1) && p[i - 1] == t[j - 1]) {
            // The last operation has been a "pass", meaning there was no editing done
            p1.push_back(p[i - 1]);
            t1.push_back(t[j - 1]);
            --i; --j;
        } else if (_table(i, j) == _table(i - 1, j - 1) + 1 && p[i - 1] != t[j - 1]) {
            // The last several (!) operations have been substitutions
            p1.push_back(']');
            t1.push_back(']');
            while (i && j && _table(i, j) == _table(i - 1, j - 1) + 1 && p[i - 1] != t[j - 1]) {
                p1.push_back(p[i - 1]);
                t1.push_back(t[j - 1]);
                --i; --j;
            }
            p1.push_back('[');
            t1.push_back('[');
        } else if (_table(i, j) == _table(i, j - 1) + 1) {
            // The last operation has been insertion into p
            p1.push_back('-');
            t1.push_back(t[j - 1]);
            --j;
        } else /*_table(i, j) == _table(i, j - 1) + 1*/ {
            // The last operation has been deletion from p
            p1.push_back(p[i - 1]);
            t1.push_back('-');
            --i;
        }
    }
    // Reaching either the first row or the first column means the rest
    // of the operations have been only insertions or only deletions.
    while (i > 0) {
        p1.push_back(p[i - 1]);
        t1.push_back('-');
        --i;
    }
    while (j > 0) {
        p1.push_back('-');
        t1.push_back(t[j - 1]);
        --j;
    }
    std::reverse(p1.begin(), p1.end());
    std::reverse(t1.begin(), t1.end());
    std::cout << "Edit process:\n" << p1 << "\n" << t1 << "\n\n";
#undef _table // for good measure
}

// Approximate string matching: the full table
void fullDynProgDemo(std::string_view p, std::string_view t) {
    const size_t m = p.length() + 1, n = t.length() + 1;
    // simulating a matrix with table[i][j] <=> array[i*n + j]
    std::vector<size_t> table(m*n);
#define _table(i,j) table[(i)*n + (j)]
    for (size_t i = 0; i < m; i++) _table(i, 0) = i;
    for (size_t j = 0; j < n; j++) _table(0, j) = 0;

    for (size_t i = 1; i < m; i++)
        for (size_t j = 1; j < n; j++) {
            _table(i, j) = min(_table(i - 1, j) + 1,
                               _table(i, j - 1) + 1,
                               _table(i - 1, j - 1) + 1 - (p[i - 1] == t[j - 1]));
        }

    // return _table(m - 1, n - 1);
    // Pretty printing the entire matrix
    const size_t numw = std::to_string(m > n ? m : n).length();
    std::cout << std::string(numw + 2, ' ');
    for (size_t j = 0; j < n - 1; j++)
        std::cout << " " << std::setw(numw) << t[j];
    std::cout << "\n";
    for (size_t i = 0; i < m; i++) {
        std::cout << ((i == 0) ? ' ' : p[i - 1]);
        for (size_t j = 0; j < n; j++)
            std::cout << " " << std::setw(numw) << _table(i, j);
        std::cout << "\n";
    }
    size_t minDist = _table(m - 1, 0);
    for (size_t j = 1; j < n; j++)
        minDist = min(minDist, _table(m - 1, j));
    std::vector<size_t> idxs;
    for (size_t j = 0; j < n; j++)
        if (_table(m - 1, j) == minDist)
            idxs.push_back(j);
    std::cout << "Minimum edit distance = " << minDist << "\n";
    std::cout << "pattern found in positions:";
    for (auto idx : idxs)
        std::cout << ' ' << idx;
    std::cout << "\n\n";
#undef _table
}

// Approximate string matching:
// The classic, space-saving alrogithm: T(m,n) = O(mn); M(m,n) = O(m)
size_t fullDynProg(std::string_view p, std::string_view t, const size_t k) {
    const size_t m = p.length(), n = t.length();
    if (k >= m)
        return 0;
    std::vector<size_t> previous(m + 1);
    std::vector<size_t> current(m + 1);
    for (size_t i = 0; i <= m; i++)
        previous[i] = i;
    size_t matches = 0;

    for (size_t i = 0; i < n; i++) {
        current[0] = 0;
        for (size_t j = 1; j <= m; j++) {
            current[j] = min(current[j - 1] + 1,
                             previous[j] + 1,
                             previous[j - 1] + 1 - (t[i] == p[j - 1]));
        }
        if (current[m] <= k)
            ++matches; // matching substring found at position i+1
        std::swap(previous, current);
    }
    return matches;
}

// Additional time optimisation - not calculating every column entirely:
// T(m,n) = O(kn) average and O(mn) worst-case; M(m,n) = O(m)
size_t cutOffDynProg(std::string_view p, std::string_view t, const size_t k) {
    const size_t m = p.length(), n = t.length();
    if (k >= m)
        return 0;
    std::vector<size_t> previous(m + 1);
    std::vector<size_t> current(m + 1);
    for (size_t i = 0; i <= k; i++)
        previous[i] = i;
    size_t matches = 0;
    size_t lp = k, lc, currSize;

    for (size_t i = 1; i <= n; i++) {
        currSize = min(lp + 1, m) + 1;
        current[0] = 0;
        lc = 0;
        for (size_t j = 1; j < currSize; j++) {
            current[j] = min(current[j - 1] + 1, previous[j - 1] + 1 - (t[i - 1] == p[j - 1]));
            if (j <= lp)
                current[j] = min(current[j], previous[j] + 1);
            if (current[j] <= k)
                lc = j;
        }
        if (currSize == m + 1 && current[m] <= k)
            ++matches; // matching substring found at position i
        std::swap(previous, current);
        lp = lc;
    }
    return matches;
}

// Calculating only certain parts from each diagonal:
// T(m,n) = O(kn) average and O(mn) worst-case; M(m,n) = O(k)
size_t diagTransitions(std::string_view p, std::string_view t, const size_t k) {
    // signed integer type needed for representation of -infinity
    using int_t = std::make_signed_t<size_t>;
    const int_t m = p.length(), n = t.length();
    if (k >= p.length())
        return 0;

    int_t table_w = n - m + size_t(k) + 3;
    int_t table_h = k + 2;
    int_t col, d;
    std::vector<int_t> preprevious(table_h, std::numeric_limits<int_t>::min());
    std::vector<int_t> previous(table_h, -1);
    std::vector<int_t> current(table_h);

    size_t matches = 0;

    for (int_t j = 2; j < table_w; j++) {
        current[0] = j - 2;
        for (int_t i = 1; i < table_h; i++) {
            col = max(preprevious[i - 1] + 1,
                      previous[i - 1] + 1,
                      current[i - 1]);
            d = j - i - 1;
            while (col < n && col - d < m && p[col - d] == t[col])
                ++col;
            current[i] = min(col, m + d);
        }
        if (current[k + 1] == m + j - k - 2)
            ++matches; // matching substring found at position current[k + 1]
        std::swap(preprevious, previous);
        std::swap(previous, current);
    }
    return matches;
}

// Bitwise parallelized calculation of an entire column. When alphabet size = b (in practise O(1)) we have:
// T(m,n) = O(bm + n); M(n) = O(bm)
size_t bitParallelized(std::string_view p, std::string_view t, const size_t k) {
    using bitmask = uint64_t;
    const size_t m = p.length(), n = t.length();
    if (k >= m || m > CHAR_BIT * sizeof(bitmask))
        return 0;
    bitmask pv = std::numeric_limits<bitmask>::max(), mv = 0;
    bitmask eq, xv, xh, ph, mh;
    size_t score = m, matches = 0;
    bitmask peq[std::numeric_limits<char>::max()] = { 0 };

    for (size_t i = 0; i < m; i++)
        peq[p[i]] |= 1ui64 << i;

    for (size_t i = 1; i <= n; i++) {
        eq = peq[t[i - 1]];
        xv = eq | mv;
        xh = (((eq & pv) + pv) ^ pv) | eq;
        ph = mv | (~(xh | pv));
        mh = pv & xh;
        if (ph & (1ui64 << (m - 1))) ++score;
        else if (mh & (1ui64 << (m - 1))) --score;
        ph <<= 1;
        pv = (mh << 1) | (~(xv | ph));
        mv = ph & xv;
        if (score <= k)
            ++matches; // matching substring found at position i
    }
    return matches;
}

int main() {
    //levenshteinDemo("polynomial", "exponential");
    //levenshteinDemo("bananas", "pyjamas");
    //// When allowing transpositions affects the result:
    //levenshteinDemo("andi", "nadq");
    //levenshteinDemo("andi", "nadq", true);

    while (true) { // Simply change this to while (false) in order to skip the interactive loop
        std::string pattern;
        std::string text;
        bool tr;
        auto input = [](auto& x) { std::cout << "> "; std::cin >> x; };
        input(pattern);
        input(text);
        input(tr);
        levenshteinDemo(pattern, text, tr);
        //fullDynProgDemo(pattern, text);
    }

    fullDynProgDemo("abc", "xyzxyacxyzxyzxyzbxyzxyz");

    // Simple sanity check:
    std::string_view p1 = "survey", p2 = "XYZsurgeryXYZXYZXYZsurgeryXYZ";
    for (auto f : { fullDynProg,cutOffDynProg,diagTransitions,bitParallelized })
        std::cout << f(p1, p2, 2) << " ";
    std::cout << "\n";
}
