/**
 *  \date April 20, 2010
 *  CTCC, University of Troms
 *
 */

#include <set>
#include <iomanip>
#include <config.h>
#include <boost/timer.hpp>

#include "constants.h"
#include "MWTree.h"

using namespace std;
using namespace Eigen;

template<int D> int MWTree<D>::defaultSplitType = NormalSplit;
template<int D> int MWTree<D>::defaultScalingType = Interpol;

/** MWTree constructor.
  * Creates an empty tree object. Node construction and assignment of most of
  * the parameters are done in derived classes. */
template<int D>
MWTree<D>::MWTree(int type, int k, const BoundingBox<D> *box)
        : MRTree<D>(k, box) {
    this->squareNorm = 0.0;
    this->autoCleanGenerated = true;

//    setupScalingBasis(type);
//    setupFilters(type);

    allocWorkMemory();
}

/** MWTree copy constructor.
  * Takes the parameters of the input tree, not it's data */
template<int D>
MWTree<D>::MWTree(const MWTree<D> &tree) : MRTree<D>(tree) {
    NOT_IMPLEMENTED_ABORT;
    this->squareNorm = 0.0;
    this->nAllocGenNodes = 0;

//    setupScalingBasis(tree.scalingType);
//    setupFilters(tree.scalingType);

    allocWorkMemory();
}

template<int D>
MWTree<D>& MWTree<D>::operator=(const MWTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

/** MWTree destructor. */
template<int D>
MWTree<D>::~MWTree() {
    this->freeWorkMemory();
}

/** Allocate work memory of the tree, for use in mwTransform and the like. */
template<int D>
void MWTree<D>::allocWorkMemory() {
    this->tmpCoefs = new Eigen::MatrixXd *[this->nThreads];
    this->tmpVector = new Eigen::VectorXd *[this->nThreads];
    this->tmpMWCoefs = new Eigen::VectorXd *[this->nThreads];
    this->nGenNodes = new int[this->nThreads];
    this->nAllocGenNodes = new int[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) {
        this->tmpCoefs[i] = new Eigen::MatrixXd(this->order + 1, D);
        this->tmpVector[i] = new Eigen::VectorXd(this->kp1_d);
        this->tmpMWCoefs[i] = new Eigen::VectorXd(this->kp1_d * (1 << D));
        this->nGenNodes[i] = 0;
        this->nAllocGenNodes[i] = 0;
    }
}

/** Deallocate work memory */
template<int D>
void MWTree<D>::freeWorkMemory() {
    for (int i = 0; i < this->nThreads; i++) {
        delete this->tmpCoefs[i];
        delete this->tmpVector[i];
        delete this->tmpMWCoefs[i];
    }
    delete[] this->tmpCoefs;
    delete[] this->tmpVector;
    delete[] this->tmpMWCoefs;
}

template<int D>
void MWTree<D>::setDefaultSplitType(int type) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void MWTree<D>::setDefaultScalingType(int type) {
    NOT_IMPLEMENTED_ABORT
}

/** Set tree split type.
  * Determines the strictness of the accuracy criterion of the wavelet norm. */
/*
template<int D>
void MWTree<D>::setSplitType(int type) {
    switch (type) {
    case (ExactSplit):
        this->splitType = ExactSplit;
        break;
    case (QuickSplit):
        this->splitType = QuickSplit;
        break;
    case (FastSplit):
        this->splitType = FastSplit;
        break;
    default:
        THROW_ERROR("Invalid tree type: " << type)
    }
}
*/

/** Calculate the squared norm of a function represented as a tree.
  *
  * Norm is calculated using endNodes only, but if your endNodeTable is
  * incomplete (e.g. within growTree), the missing nodes must be given in the
  * input work vector. Involves an MPI reduction operation. */
template<int D>
double MWTree<D>::calcTreeNorm(MRNodeVector *work)  {
    NOT_IMPLEMENTED_ABORT;
//    double treeNorm = 0.0;
//    int nNodes = 0;
//    if (work != 0) {
//        nNodes = work->size();
//    }
//    for (int n = 0; n < nNodes; n++) {
//        MWNode<D> *node = (*work)[n];
//        if (not node->isForeign()) {
//            assert(node->hasCoefs());
//            treeNorm += node->getSquareNorm();
//        }
//    }
//    for (int n = 0; n < (int) this->endNodeTable.size(); n++) {
//        MWNode<D> *node = this->endNodeTable[n];
//        if ( not node->isForeign()) {
//            treeNorm += node->getSquareNorm();
//        }
//    }
//    treeNorm = reduceNorm(treeNorm);
//    return treeNorm;
}

/** Collect local treeNorms to global treeNorm. */
template<int D>
double MWTree<D>::reduceNorm(double treeNorm) {
    NOT_IMPLEMENTED_ABORT;
//#ifdef HAVE_MPI
//    if (isBuildDistributed()) {
//        return mpi::all_reduce(node_group, treeNorm, std::plus<double>());
//    }
//#endif
//    return treeNorm;
}

/** Reduce the accuracy of the tree by deleting nodes
  * which have a higher precision than the requested precison.
  * By default, the relative precision of the tree is used. */
template<int D>
void MWTree<D>::cropTree(double prec, bool absPrec) {
    NOT_IMPLEMENTED_ABORT;
//    set<const NodeIndex<D> *, NodeIndexComp<D> > cropNodes;
//    for (int i = 0; i < this->rootBox.getNBoxes(); i++) {
//        MWNode<D> &rootNode = getRootMWNode(i);
//        if (this->isScattered()) {
//            rootNode.cropChildren(prec, &cropNodes);
//        } else {
//            rootNode.cropChildren(prec);
//        }
//    }
//    if (this->isScattered()) {
//        broadcastIndexList(cropNodes);
//        typename set<const NodeIndex<D> *,
//                NodeIndexComp<D> >::reverse_iterator it;
//        for (it = cropNodes.rbegin(); it != cropNodes.rend(); it++) {
//            MWNode<D> *node = this->findNode(**it);
//            if (node != 0) {
//                node->deleteChildren();
//            }
//        }
//    }
//    resetEndNodeTable();
//    this->squareNorm = calcTreeNorm();
}

/** Regenerate all s/d-coeffs by backtransformation, starting at the bottom and
  * thus purifying all coefficients. Option to overwrite or add up existing
  * coefficients of BranchNodes (can be used after operator application). */
template<int D>
void MWTree<D>::mwTransformUp(bool overwrite) {
    NOT_IMPLEMENTED_ABORT;
//    vector<MWNodeVector > nodeTable;
//    if (isScattered()) {
//        this->makeLocalNodeTable(nodeTable);
//    } else {
//        makeNodeTable(nodeTable);
//    }
//    set<MWNode<D> *> missing;
//    int start = nodeTable.size() - 2;
//    for (int n = start; n >= 0; n--) {
//        if (isScattered()) {
//            missing.clear();
//            findMissingChildren(nodeTable[n], missing);
//            //communicate missing
//            syncNodes(missing);
//        }
//        int n_nodes = nodeTable[n].size();
//#pragma omp parallel firstprivate(n_nodes, overwrite) shared(nodeTable)
//        {
//#pragma omp for schedule(guided)
//            for (int i = 0; i < n_nodes; i++) {
//                MWNode<D> &node = *nodeTable[n][i];
//                if (node.isBranchNode()) {
//                    node.reCompress(overwrite);
//                }
//            }
//        }
//        // For some very strange reason the calcNorms routine cannot be called in the
//        // above OMP loop.
//        for (int i = 0; i < n_nodes; i++) {
//            MWNode<D> &node = *nodeTable[n][i];
//            if (node.isBranchNode()) {
//                node.calcNorms();
//            }
//        }
//    }
}

/** Regenerate all scaling coeffs by MW transformation of existing s/w-coeffs
  * on coarser scales, starting at the rootNodes. Option to overwrite or add up
  * existing scaling coefficients (can be used after operator application). */
template<int D>
void MWTree<D>::mwTransformDown(bool overwrite) {
    NOT_IMPLEMENTED_ABORT;
//    this->purgeForeignNodes();
//    if (isScattered()) {
//        vector<MWNodeVector > nodeTable;
//        this->makeLocalNodeTable(nodeTable);

//        vector<MWNodeVector > parentTable;
//        this->makeNodeTable(parentTable);

//        set<MWNode<D> *> missing;
//        for (int i = 1; i < nodeTable.size(); i++) {
//            missing.clear();
//            findMissingParents(nodeTable[i], missing);
//            //communicate missing
//            syncNodes(missing);
//            for (unsigned int j = 0; j < parentTable[i-1].size(); j++) {
//                MWNode<D> &parent = *parentTable[i-1][j];
//                if (not parent.hasCoefs() or parent.isLeafNode()) {
//                    continue;
//                }
//                parent.giveChildrenScaling(overwrite);
//                for (int k = 0; k < parent.getNChildren(); k++) {
//                    MWNode<D> &child = parent.getMWChild(k);
//                    if (child.isForeign()) {
//                        child.clearCoefs();
//                        child.clearNorms();
//                    }
//                }
//            }
//        }
//    } else {
//        vector<MWNodeVector > nodeTable;
//        makeNodeTable(nodeTable);
//#pragma omp parallel shared(nodeTable)
//        {
//            for (int n = 0; n < nodeTable.size(); n++) {
//                int n_nodes = nodeTable[n].size();
//#pragma omp for schedule(guided)
//                for (int i = 0; i < n_nodes; i++) {
//                    MWNode<D> &node = *nodeTable[n][i];
//                    if (node.isBranchNode()) {
//                        node.giveChildrenScaling(overwrite);
//                    }
//                }
//            }
//        }
//        this->purgeForeignNodes();
//    }
}

/** Initialize the work table with the root nodes (for initial
  * projection without seeding) */
/*
template<int D>
void MWTree<D>::setupWorkTable(MWNodeVector &wt) {
    int nbox = this->rootBox.getNBoxes();
    for (int i = 0; i < nbox; i++) {
        MWNode<D> &nd = this->rootBox.getNode(i);
        wt.push_back(&nd);
    }
}
*/

/** Search the nodeTable for OWN node that does not have coefficients,
  * communicate their parents and give then scaling coeffs by MW transform. */
/*
template<int D>
void MWTree<D>::updateMissingScalingPart(const MWNodeVector &nodeTable) {
    set<MWNode<D> *> missing;

    for (unsigned int i = 0; i < nodeTable.size(); i++) {
        MWNode<D> &node = *nodeTable[i];
        if (not node.isForeign() and (not node.hasCoefs())) {
            missing.insert(&node.getMWParent());
        }
    }
    syncNodes(missing);
    typename set<MWNode<D> *>::iterator it;
    for (it = missing.begin(); it != missing.end(); ++it) {
        (*it)->giveChildrenScaling();
    }
}
*/

/** Tag each node with the rank who owns it. */
/*
template<int D>
void MWTree<D>::assignNodeTags(MWNodeVector &workTable, int rank) {

    int nNodes = workTable.size();
    int start = 0;
    int end = nNodes;
    int nLocales = 1;

    if (isBuildDistributed()) {
        nLocales = node_group.size();
    }
*/
    /* The default rank is -1, which tags nodes according to which locale
        they belong to, if the tree is distributed */
/*
    int loc = rank;
    if (rank < 0) {
        loc = this->getRankId();
    }
    for (int k = 0; k < nLocales; k++) {
        if (isBuildDistributed()) {
            if (rank < 0) {
                loc = k;
            }
            get_locale_index_range(k, nNodes, start, end);
        }

        for (int i = start; i < end; i++) {
            MWNode<D> &node = *workTable[i];
            node.setRankId(loc);
        }
    }
}
*/

/** Initialize MW filter cache. */
template<int D>
void MWTree<D>::setupFilters(int type) {
    NOT_IMPLEMENTED_ABORT;
//    getLegendreFilterCache(lfilters);
//    getInterpolatingFilterCache(ifilters);
//    switch (type) {
//    case Legendre:
//        this->filter = &lfilters.get(this->order);
//        break;
//    case Interpol:
//        this->filter = &ifilters.get(this->order);
//        break;
//    default:
//        THROW_ERROR("Invalid scaling basis selected.")
//    }
}

/** Initialize scaling basis cache. */
template<int D>
void MWTree<D>::setupScalingBasis(int type) {
    NOT_IMPLEMENTED_ABORT;
//    this->scalingType = type;
//    getLegendreScalingCache(lsf);
//    getInterpolatingScalingCache(isf);
//    switch (type) {
//    case Legendre:
//        this->scalingFunc = &lsf.get(this->order);
//        break;
//    case Interpol:
//        this->scalingFunc = &isf.get(this->order);
//        break;
//    default:
//        THROW_ERROR("Invalid scaling basis selected.")
//    }
}

/** Traverse tree and set all nodes to zero.
  *
  * Keeps the node structure of the tree, even though the zero function is
  * representable at depth zero. Use cropTree to remove unnecessary nodes.
  * Option to keep treeNorm (used after seed). */
template<int D>
void MWTree<D>::setZero(bool clearTreeNorm) {
    NOT_IMPLEMENTED_ABORT;
//    purgeForeignNodes(false);
//    LebesgueIterator<D> it(this);
//    while(it.next()) {
//        MWNode<D> &node = it.getNode();
//        if (not node.isForeign()) {
//            node.getCoefs().setZero();
//            node.calcNorms();
//        }
//    }
//    if (clearTreeNorm) {
//        this->squareNorm = 0.0;
//    }
}


/** Traverse tree and set all nodeWeights to zero. */
//template<int D>
//void MWTree<D>::clearNodeWeights() {
//    HilbertIterator<D> it(this);
//    while (it.next()) {
//        MWNode<D> &node = it.getNode();
//        node.clearNodeWeight(0);
//        node.clearNodeWeight(1);
//    }
//}

template class MWTree<1>;
template class MWTree<2>;
template class MWTree<3>;
