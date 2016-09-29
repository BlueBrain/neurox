#include <stdio.h>
#include "coreneuron/nrniv/nrn_assert.h"
#include "coreneuron/nrniv/cellorder.h"
#include "coreneuron/nrniv/tnode.h"

// just for use_interleave_permute
#include "coreneuron/nrniv/nrniv_decl.h"

#include <map>
#include <set>
#include <algorithm>
#include <string.h>

using namespace std;

static size_t groupsize = 32;

static bool tnode_earlier(TNode* a, TNode* b) {
    bool result = false;
    if (a->treesize < b->treesize) {  // treesize dominates
        result = true;
    } else if (a->treesize == b->treesize) {
        if (a->hash < b->hash) {  // if treesize same, keep identical trees together
            result = true;
        } else if (a->hash == b->hash) {
            result = a->nodeindex < b->nodeindex;  // identical trees ordered by nodeindex
        }
    }
    return result;
}

static bool ptr_tnode_earlier(TNode* a, TNode* b) {
    return tnode_earlier(a, b);
}

TNode::TNode(int ix) {
    nodeindex = ix;
    cellindex = 0;
    groupindex = 0;
    level = 0;
    hash = 0;
    treesize = 1;
    nodevec_index = 0;
    treenode_order = 0;
    parent = NULL;
    children.reserve(2);
}

TNode::~TNode() {
}

size_t TNode::mkhash() {  // call on all nodes in leaf to root order
    // concept from http://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
    std::sort(children.begin(), children.end(), ptr_tnode_earlier);
    hash = children.size();
    treesize = 1;
    for (size_t i = 0; i < children.size(); ++i) {  // need sorted by child hash
        hash ^= children[i]->hash + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        treesize += children[i]->treesize;
    }
    return hash;  // hash of leaf nodes is 0
}

static void tree_analysis(int* parent, int nnode, int ncell, VecTNode&);
static void node_interleave_order(int ncell, VecTNode&);
static void admin1(int ncell,
                   VecTNode& nodevec,
                   int& nwarp,
                   int& nstride,
                   int*& stride,
                   int*& firstnode,
                   int*& lastnode,
                   int*& cellsize);
static void admin2(int ncell,
                   VecTNode& nodevec,
                   int& nwarp,
                   int& nstride,
                   int*& stridedispl,
                   int*& strides,
                   int*& rootbegin,
                   int*& nodebegin,
                   int*& ncycles);
static void check(VecTNode&);
static void prtree(VecTNode&);

typedef std::pair<TNode*, int> TNI;
typedef std::map<size_t, pair<TNode*, int> > HashCnt;
typedef vector<TNI> TNIVec;

static char* stree(TNode* nd) {
    char s[1000];

    if (nd->treesize > 100) {
        return strdup("");
    }
    s[0] = '(';
    s[1] = '\0';
    for (size_t i = 0; i < nd->children.size(); ++i) {  // need sorted by child hash
        char* sr = stree(nd->children[i]);
        strcat(s, sr);
        free(sr);
    }
    strcat(s, ")");
    return strdup(s);
}

/*
assess the quality of the ordering. The measure is the size of a contiguous
list of nodes whose parents have the same order. How many contiguous lists
have that same size. How many nodes participate in that size list.
Modify the quality measure from experience with performance. Start with
list of (nnode, size_participation)
*/
static void quality(VecTNode& nodevec, size_t max = 32) {
    size_t qcnt = 0;  // how many contiguous nodes have contiguous parents

    // first ncell nodes are by definition in contiguous order
    for (size_t i = 0; i < nodevec.size(); ++i) {
        if (nodevec[i]->parent != NULL) {
            break;
        }
        qcnt += 1;
    }
    size_t ncell = qcnt;

    // key is how many parents in contiguous order
    // value is number of nodes that participate in that
    map<size_t, size_t> qual;
    size_t ip_last = 10000000000;
    for (size_t i = ncell; i < nodevec.size(); ++i) {
        size_t ip = nodevec[i]->parent->nodevec_index;
        // i%max == 0 means that if we start a warp with 8 and then have 32
        // the 32 is broken into 24 and 8. (modify if the arrangement during
        // gaussian elimination becomes more sophisticated.(
        if (ip == ip_last + 1 && i % max != 0) {  // contiguous
            qcnt += 1;
        } else {
            if (qcnt == 1) {
                // printf("unique %ld p=%ld ix=%d\n", i, ip, nodevec[i]->nodeindex);
            }
            qual[max] += (qcnt / max) * max;
            size_t x = qcnt % max;
            if (x) {
                qual[x] += x;
            }
            qcnt = 1;
        }
        ip_last = ip;
    }
    qual[max] += (qcnt / max) * max;
    size_t x = qcnt % max;
    if (x) {
        qual[x] += x;
    }

    // print result
    qcnt = 0;
#if 0
  for (map<size_t, size_t>::iterator it = qual.begin(); it != qual.end(); ++it) {
    qcnt += it->second;
    printf("%6ld %6ld\n", it->first, it->second);
  }
#endif
#if 0
  printf("qual.size=%ld  qual total nodes=%ld  nodevec.size=%ld\n",
    qual.size(), qcnt, nodevec.size());
#endif

    // how many race conditions. ie refer to same parent on different core
    // of warp (max cores) or parent in same group of max.
    size_t maxip = ncell;
    size_t nrace1 = 0;
    size_t nrace2 = 0;
    set<size_t> ipused;
    for (size_t i = ncell; i < nodevec.size(); ++i) {
        TNode* nd = nodevec[i];
        size_t ip = nd->parent->nodevec_index;
        if (i % max == 0) {
            maxip = i;
            ipused.clear();
        }
        if (ip >= maxip) {
            nrace1 += 1;
        } /*else*/
        {
            if (ipused.find(ip) != ipused.end()) {
                nrace2 += 1;
                if (ip >= maxip) {
                    // printf("race for parent %ld (parent in same group as multiple users))\n",
                    // ip);
                }
            } else {
                ipused.insert(ip);
            }
        }
    }
#if 0
  printf("nrace = %ld (parent in same group of %ld nodes)\n", nrace1, max);
  printf("nrace = %ld (parent used more than once by same group of %ld nodes)\n", nrace2, max);
#endif
}

size_t level_from_root(VecTNode& nodevec) {
    size_t maxlevel = 0;
    for (size_t i = 0; i < nodevec.size(); ++i) {
        TNode* nd = nodevec[i];
        if (nd->parent) {
            nd->level = nd->parent->level + 1;
            if (maxlevel < nd->level) {
                maxlevel = nd->level;
            }
        } else {
            nd->level = 0;
        }
    }
    return maxlevel;
}

size_t level_from_leaf(VecTNode& nodevec) {
    size_t maxlevel = 0;
    for (size_t i = nodevec.size() - 1; true; --i) {
        TNode* nd = nodevec[i];
        size_t lmax = 0;
        for (size_t ichild = 0; ichild < nd->children.size(); ++ichild) {
            if (lmax <= nd->children[ichild]->level) {
                lmax = nd->children[ichild]->level + 1;
            }
        }
        nd->level = lmax;
        if (maxlevel < lmax) {
            maxlevel = lmax;
        }
        if (i == 0) {
            break;
        }
    }
    return maxlevel;
}

static void set_cellindex(int ncell, VecTNode& nodevec) {
    for (int i = 0; i < ncell; ++i) {
        nodevec[i]->cellindex = i;
    }
    for (size_t i = 0; i < nodevec.size(); ++i) {
        TNode& nd = *nodevec[i];
        for (size_t j = 0; j < nd.children.size(); ++j) {
            TNode* cnode = nd.children[j];
            cnode->cellindex = nd.cellindex;
        }
    }
}

static void set_groupindex(VecTNode& nodevec) {
    for (size_t i = 0; i < nodevec.size(); ++i) {
        TNode* nd = nodevec[i];
        if (nd->parent) {
            nd->groupindex = nd->parent->groupindex;
        } else {
            nd->groupindex = i / groupsize;
        }
    }
}

#if 0
#define MSS MSS_ident_stat
typedef map<size_t, size_t> MSS;
static bool vsmss_comp(const pair<size_t, MSS*>& a, const pair<size_t, MSS*>& b) {
  bool result = false;
  const MSS::iterator& aa = a.second->begin();
  const MSS::iterator& bb = b.second->begin();
  if (aa->first < bb->first) {
    result = true;
  }else if (aa->first == bb->first) {
    if (aa->second < bb->second) {
      result = true;
    }
  }
  return result;
}
#endif

// how many identical trees and their levels
// print when more than one instance of a type
// reverse the sense of levels (all leaves are level 0) to get a good
// idea of the depth of identical subtrees.
static void ident_statistic(VecTNode& nodevec, size_t ncell) {
    // reverse sense of levels
    //  size_t maxlevel = level_from_leaf(nodevec);
    size_t maxlevel = level_from_root(nodevec);

    // # in each level
    vector<vector<size_t> > n_in_level(maxlevel + 1);
    for (size_t i = 0; i <= maxlevel; ++i) {
        n_in_level[i].resize(ncell / groupsize);
    }
    for (size_t i = 0; i < nodevec.size(); ++i) {
        n_in_level[nodevec[i]->level][nodevec[i]->groupindex]++;
    }
    printf("n_in_level.size = %ld\n", n_in_level.size());
    for (size_t i = 0; i < n_in_level.size(); ++i) {
        printf("%5ld\n", i);
        for (size_t j = 0; j < n_in_level[i].size(); ++j) {
            printf(" %5ld", n_in_level[i][j]);
        }
        printf("\n");
    }

#if 0
  typedef map<size_t, MSS> MSMSS;
  typedef vector<pair<size_t, MSS*> > VSMSS;
  MSMSS info;
  for (size_t i=0; i < nodevec.size(); ++i) {
    TNode* nd = nodevec[i];
    info[nd->hash][nd->level]++;
  }

  VSMSS vinfo;
  for (MSMSS::iterator i = info.begin(); i != info.end(); ++i) {
    vinfo.push_back(pair<size_t, MSS*>(i->first, &(i->second)));
  }
  std::sort(vinfo.begin(), vinfo.end(), vsmss_comp);

  for (VSMSS::iterator i = vinfo.begin(); i < vinfo.end(); ++i) {
    MSS* ival = i->second;
    if (ival->size() > 1 || ival->begin()->second > 8) {
      printf("hash %ld", i->first);
      for (MSS::iterator j = ival->begin(); j != ival->end(); ++j) {
        printf(" (%ld, %ld)", j->first, j->second);
      }
      printf("\n");
    }
  }
  printf("max level = %ld\n", maxlevel);
#endif
}
#undef MSS

// for cells with same size, keep identical trees together

// parent is (unpermuted)  nnode length vector of parent node indices.
// return a permutation (of length nnode) which orders cells of same
// size so that identical trees are grouped together.
// Note: cellorder[ncell:nnode] are the identify permutation.

int* node_order(int ncell,
                int nnode,
                int* parent,
                int& nwarp,
                int& nstride,
                int*& stride,
                int*& firstnode,
                int*& lastnode,
                int*& cellsize,
                int*& stridedispl) {
    VecTNode nodevec;
    if (0)
        prtree(nodevec);  // avoid unused warning

    // nodevec[0:ncell] in increasing size, with identical trees together,
    // and otherwise nodeindex order
    tree_analysis(parent, nnode, ncell, nodevec);
    check(nodevec);

    set_cellindex(ncell, nodevec);
    set_groupindex(nodevec);
    level_from_root(nodevec);

    // nodevec[ncell:nnode] cells are interleaved in nodevec[0:ncell] cell order
    if (use_interleave_permute == 1) {
        node_interleave_order(ncell, nodevec);
    } else {
        group_order2(nodevec, groupsize, ncell);
    }
    check(nodevec);

#if 0
  for (int i=0; i < ncell; ++i) {
    TNode& nd = *nodevec[i];
    printf("%d size=%ld hash=%ld ix=%d\n", i, nd.treesize, nd.hash, nd.nodeindex);
  }
#endif

    if (0)
        ident_statistic(nodevec, ncell);
    quality(nodevec);

    // the permutation
    int* nodeorder = new int[nnode];
    for (int i = 0; i < nnode; ++i) {
        TNode& nd = *nodevec[i];
        nodeorder[nd.nodeindex] = i;
    }

    // administrative statistics for gauss elimination
    if (use_interleave_permute == 1) {
        admin1(ncell, nodevec, nwarp, nstride, stride, firstnode, lastnode, cellsize);
    } else {
        //  admin2(ncell, nodevec, nwarp, nstride, stridedispl, stride, rootbegin, nodebegin,
        //  ncycles);
        admin2(ncell, nodevec, nwarp, nstride, stridedispl, stride, firstnode, lastnode, cellsize);
    }

#if 1
    int ntopol = 1;
    for (int i = 1; i < ncell; ++i) {
        if (nodevec[i - 1]->hash != nodevec[i]->hash) {
            ntopol += 1;
        }
    }
    printf("%d distinct tree topologies\n", ntopol);
#endif

    for (size_t i = 0; i < nodevec.size(); ++i) {
        delete nodevec[i];
    }

    return nodeorder;
}

void check(VecTNode& nodevec) {
    // printf("check\n");
    size_t nnode = nodevec.size();
    size_t ncell = 0;
    for (size_t i = 0; i < nnode; ++i) {
        nodevec[i]->nodevec_index = i;
        if (nodevec[i]->parent == NULL) {
            ncell++;
        }
    }
    for (size_t i = 0; i < ncell; ++i) {
        nrn_assert(nodevec[i]->parent == NULL);
    }
    for (size_t i = ncell; i < nnode; ++i) {
        TNode& nd = *nodevec[i];
        if (nd.parent->nodevec_index >= nd.nodevec_index) {
            printf("error i=%ld nodevec_index=%ld parent=%ld\n", i, nd.nodevec_index,
                   nd.parent->nodevec_index);
        }
        nrn_assert(nd.nodevec_index > nd.parent->nodevec_index);
    }
}

void prtree(VecTNode& nodevec) {
    size_t nnode = nodevec.size();
    for (size_t i = 0; i < nnode; ++i) {
        nodevec[i]->nodevec_index = i;
    }
    for (size_t i = 0; i < nnode; ++i) {
        TNode& nd = *nodevec[i];
        printf("%ld p=%d   c=%ld l=%ld o=%ld   ix=%d pix=%d\n", i,
               nd.parent ? int(nd.parent->nodevec_index) : -1, nd.cellindex, nd.level,
               nd.treenode_order, nd.nodeindex, nd.parent ? int(nd.parent->nodeindex) : -1);
    }
}

void tree_analysis(int* parent, int nnode, int ncell, VecTNode& nodevec) {
    //  VecTNode nodevec;

    // create empty TNodes (knowing only their index)
    nodevec.reserve(nnode);
    for (int i = 0; i < nnode; ++i) {
        nodevec.push_back(new TNode(i));
    }

    // determine the (sorted by hash) children of each node
    for (int i = nnode - 1; i >= ncell; --i) {
        nodevec[i]->parent = nodevec[parent[i]];
        nodevec[i]->mkhash();
        nodevec[parent[i]]->children.push_back(nodevec[i]);
    }

    // determine hash of the cells
    for (int i = 0; i < ncell; ++i) {
        nodevec[i]->mkhash();
    }

    std::sort(nodevec.begin(), nodevec.begin() + ncell, tnode_earlier);
}

static bool interleave_comp(TNode* a, TNode* b) {
    bool result = false;
    if (a->treenode_order < b->treenode_order) {
        result = true;
    } else if (a->treenode_order == b->treenode_order) {
        if (a->cellindex < b->cellindex) {
            result = true;
        }
    }
    return result;
}

// sort so nodevec[ncell:nnode] cell instances are interleaved. Keep the
// secondary ordering with respect to treenode_order so each cell is still a tree.

void node_interleave_order(int ncell, VecTNode& nodevec) {
    int* order = new int[ncell];
    for (int i = 0; i < ncell; ++i) {
        order[i] = 0;
        nodevec[i]->treenode_order = order[i]++;
    }
    for (size_t i = 0; i < nodevec.size(); ++i) {
        TNode& nd = *nodevec[i];
        for (size_t j = 0; j < nd.children.size(); ++j) {
            TNode* cnode = nd.children[j];
            cnode->treenode_order = order[nd.cellindex]++;
        }
    }
    delete[] order;

    //  std::sort(nodevec.begin() + ncell, nodevec.end(), contig_comp);
    std::sort(nodevec.begin() + ncell, nodevec.end(), interleave_comp);

#if 0
  for (size_t i=0; i < nodevec.size(); ++i) {
    TNode& nd = *nodevec[i];
    printf("%ld cell=%ld ix=%d\n",  i, nd.cellindex, nd.nodeindex);
  }
#endif
}

static void admin1(int ncell,
                   VecTNode& nodevec,
                   int& nwarp,
                   int& nstride,
                   int*& stride,
                   int*& firstnode,
                   int*& lastnode,
                   int*& cellsize) {
    // firstnode[i] is the index of the first nonroot node of the cell
    // lastnode[i] is the index of the last node of the cell
    // cellsize is the number of nodes in the cell not counting root.
    // nstride is the maximum cell size (not counting root)
    // stride[i] is the number of cells with an ith node.
    firstnode = new int[ncell];
    lastnode = new int[ncell];
    cellsize = new int[ncell];

    nwarp = (ncell % warpsize == 0) ? (ncell / warpsize) : (ncell / warpsize + 1);

    for (int i = 0; i < ncell; ++i) {
        firstnode[i] = -1;
        lastnode[i] = -1;
        cellsize[i] = 0;
    }

    nstride = 0;
    for (size_t i = ncell; i < nodevec.size(); ++i) {
        TNode& nd = *nodevec[i];
        size_t ci = nd.cellindex;
        if (firstnode[ci] == -1) {
            firstnode[ci] = i;
        }
        lastnode[ci] = i;
        cellsize[ci] += 1;
        if (nstride < cellsize[ci]) {
            nstride = cellsize[ci];
        }
    }

    stride = new int[nstride + 1];  // in case back substitution accesses this
    for (int i = 0; i <= nstride; ++i) {
        stride[i] = 0;
    }
    for (size_t i = ncell; i < nodevec.size(); ++i) {
        TNode& nd = *nodevec[i];
        stride[nd.treenode_order - 1] += 1;  // -1 because treenode order includes root
    }
}

// for admin2 we allow the node organisation in warps of (say 4 cores per warp)
// ...............  ideal warp but unbalanced relative to warp with max cycles
// ...............  ncycle = 15, icore [0:4), all strides are 4.
// ...............
// ...............
//
// ..........       unbalanced relative to warp with max cycles
// ..........       ncycle = 10, not all strides the same because
// ..........       of need to avoid occasional race conditions.
//  .  . ..         icore [4:8) only 4 strides of 4
//
// ....................  ncycle = 20, uses only one core in the warp (cable)
//                       icore 8, all ncycle strides are 1

// One thing to be unhappy about is the large stride vector of size about
// number of compartments/warpsize. There are a lot of models where the
// stride for a warp is constant except for one cycle in the warp and that
// is easy to obtain when there are more than warpsize cells per warp.

static size_t stride_length(size_t begin, size_t end, VecTNode& nodevec) {
    // return stride length starting at i. Do not go past j.
    // max stride is warpsize.
    // At this time, only assume vicious parent race conditions matter.
    if (end - begin > warpsize) {
        end = begin + warpsize;
    }
    for (size_t i = begin; i < end; ++i) {
        TNode* nd = nodevec[i];
        nrn_assert(nd->nodevec_index == i);
        size_t diff = dist2child(nd);
        if (i + diff < end) {
            end = i + diff;
        }
    }
    return end - begin;
}

static void admin2(int ncell,
                   VecTNode& nodevec,
                   int& nwarp,
                   int& nstride,
                   int*& stridedispl,
                   int*& strides,
                   int*& rootbegin,
                   int*& nodebegin,
                   int*& ncycles) {
    // the number of groups is the number of warps needed
    // ncore is the number of warps * warpsize
    nwarp = nodevec[ncell - 1]->groupindex + 1;

    ncycles = new int[nwarp];
    stridedispl = new int[nwarp + 1];  // running sum of ncycles (start at 0)
    rootbegin = new int[nwarp + 1];    // index (+1) of first root in warp.
    nodebegin = new int[nwarp + 1];    // index (+1) of first node in warp.

    // rootbegin and nodebegin are the root index values + 1 of the last of
    // the sequence of constant groupindex
    rootbegin[0] = 0;
    for (size_t i = 0; i < size_t(ncell); ++i) {
        rootbegin[nodevec[i]->groupindex + 1] = i + 1;
    }
    nodebegin[0] = ncell;
    for (size_t i = size_t(ncell); i < nodevec.size(); ++i) {
        nodebegin[nodevec[i]->groupindex + 1] = i + 1;
    }

    // ncycles, stridedispl, and nstride
    nstride = 0;
    stridedispl[0] = 0;
    for (size_t iwarp = 0; iwarp < (size_t)nwarp; ++iwarp) {
        size_t j = size_t(nodebegin[iwarp + 1]);
        int nc = 0;
        size_t i = nodebegin[iwarp];
        while (i < j) {
            i += stride_length(i, j, nodevec);
            ++nc;
        }
        ncycles[iwarp] = nc;
        stridedispl[iwarp + 1] = stridedispl[iwarp] + nc;
        nstride += nc;
    }

    // strides
    strides = new int[nstride];
    nstride = 0;
    for (size_t iwarp = 0; iwarp < (size_t)nwarp; ++iwarp) {
        size_t j = size_t(nodebegin[iwarp + 1]);
        size_t i = nodebegin[iwarp];
        while (i < j) {
            int k = stride_length(i, j, nodevec);
            i += k;
            strides[nstride++] = k;
        }
    }

#if 0
printf("warp rootbegin nodebegin stridedispl\n");
for (int i = 0; i <= nwarp; ++i){
  printf("%4d %4d %4d %4d\n", i, rootbegin[i], nodebegin[i], stridedispl[i]);
}
#endif
}
