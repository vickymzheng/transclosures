/***
 *  $Id$
 **
 *  File: nano-tools-stream.cpp
 *  Created: Mar 16, 2018
 *
 *  Author: Double Blind 
 *  Copyright (c) 2018 Double Blind
 *  Distributed under the [LICENSE].
 *  See accompanying file LICENSE.
 *
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <vector>

#include "parallel-hashmap/parallel_hashmap/phmap.h"


struct edge {
    int source;
    int target;
    int overlap;
    int slength;
    int tlength;
    char sorient;
    char torient;
}; // struct edge

inline bool operator==(const edge& lhs, const edge& rhs) {
    return ((lhs.source == rhs.source) && (lhs.target == rhs.target));
} // operator==

inline std::istream& operator>>(std::istream& is, edge& e) {
    is >> e.source;

    if (e.source == -1) {
        e.target = -1;
        return is;
    }

    is >> e.target >> e.overlap >> e.slength >> e.tlength;
    std::string s;
    is >> s;
    e.sorient = s[0];
    is >> s;
    e.torient = s[0];
    return is;
} // operator>>

inline std::ostream& operator<<(std::ostream& os, const edge& e) {
    os << e.source << " " << e.target << " " << e.overlap
       << " " << e.slength << " " << e.tlength
       << " " << e.sorient << " " << e.torient;
    return os;
} // operator<<


struct edge_hash {
    std::size_t operator()(const edge& e) const {
        int s = e.source;
        int t = e.target;
        return (((s + t) * (s + t + 1)) >> 1) + t;
    }
}; // struct edge_hash

edge make_edge(int s, int t) {
    edge e;
    if (t < s) std::swap(s, t);
    e.source = s;
    e.target = t;
    return e;
} // make_edge


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "usage: " << argv[0] << " edge_file seed" << std::endl;
        std::cout << "       -1: random seed, 0: no randomization" << std::endl;
        return -1;
    }

    int seed = std::stoi(argv[2]);

    // currently we read only edges
    // this can be extended later
    int n = 0;

    std::vector<std::vector<int>> G;
    phmap::flat_hash_set<edge, edge_hash> edges;

    std::ifstream f(argv[1]);

    if (!f) {
        std::cout << "error: could not read " << argv[1] << std::endl;
        return -1;
    }

    std::istream_iterator<edge> f_it(f);
    std::istream_iterator<edge> f_end;

    for (; f_it != f_end; ++f_it) {
        edge e = *f_it;
        if (e.source == -1) continue;

        n = std::max({n, e.source, e.target});
        if (G.size() < n + 1) G.resize(n + 1);

        G[e.source].push_back(e.target);
        G[e.target].push_back(e.source);

        if (e.source > e.target) std::swap(e.source, e.target);

        edges.insert(e);
    }

    f.close();

    // now we can simulate
    // first some random order of incoming reads
    n = n + 1;
    std::vector<int> ev(n);

    std::random_device rd;
    std::mt19937 gen((seed == -1) ? rd() : seed);

    std::iota(std::begin(ev), std::end(ev), 0);
    if (seed != 0) std::shuffle(std::begin(ev), std::end(ev), gen);

    // here go requests
    std::vector<bool> ev_map(n, false);

    for (auto s : ev) {
        ev_map[s] = true;

        for (auto t : G[s]) {
            if (ev_map[t]) {
                auto e_it = edges.find(make_edge(s, t));
                std::cout << *e_it << std::endl << std::flush;
            }
        } // for t

        std::cout << -1 << " " << s << std::endl << std::flush;
    } // for s

    return 0;
} // main
