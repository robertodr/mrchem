/**
 *
 *
 *  \date Aug 14, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 *
 */

#include "ProjectedNode.h"
#include "GenNode.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using namespace std;
using namespace Eigen;

/** Root node constructor. By default root nodes are initialized to
  * represent functions which are constant zero. */
template<int D>
ProjectedNode<D>::ProjectedNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx)
        : FunctionNode<D> (t, nIdx) {
    this->setIsEndNode();
    if (not this->isForeign()) {
        this->allocCoefs();
        this->zeroCoefs();
        this->zeroNorms();
    }
}

/** ProjectedNode constructor.
  * Creates an empty node given its parent and translation vector */
template<int D>
ProjectedNode<D>::ProjectedNode(ProjectedNode<D> &p, int cIdx)
        : FunctionNode<D> (p, cIdx) {
    NOT_IMPLEMENTED_ABORT;
//    this->setIsEndNode();
//    if (not this->isForeign()) {
//        this->allocCoefs();
//        this->zeroCoefs();
//        this->zeroNorms();
//    }
}

template<int D>
ProjectedNode<D>::ProjectedNode(const MWNode<D> &n)
        : FunctionNode<D>(n) {
    if (not this->isForeign()) {
        this->allocCoefs();
        this->zeroCoefs();
    }
    const VectorXd &c = n.getCoefs();
    int cSize = c.size();
    assert (cSize <= this->getNCoefs());
    this->coefs->segment(0, cSize) = c;
    this->setHasCoefs();
    this->calcNorms();
}


template<int D>
ProjectedNode<D>::ProjectedNode(const ProjectedNode<D> &n)
        : FunctionNode<D>(n) {
    NOT_IMPLEMENTED_ABORT;
//    if (not this->isForeign()) {
//        this->allocCoefs();
//        this->zeroCoefs();
//        this->zeroNorms();
//    }
}

/* Recurcive node constructor*/
template<int D>
void ProjectedNode<D>::copyChildren(const MRNode<D> &node) {
    NOT_IMPLEMENTED_ABORT;
//    if (node.isBranchNode()) {
//        this->allocKindergarten();
//        this->setIsBranchNode();
//        this->clearIsEndNode();
//    }
//    int myRank = this->getRankId();
//    for (int cIdx = 0; cIdx < node.getNChildren(); cIdx++) {
//        const MRNode<D> &yourChild = node.getMRChild(cIdx);
//        int childRank = yourChild.getRankId();
//        this->setRankId(childRank); //Rank is copied from parent
//        ProjectedNode<D> *myChild = new ProjectedNode(*this, cIdx);
//        this->setRankId(myRank);
//        myChild->copyChildren(yourChild);
//        this->children[cIdx] = myChild;
//    }
}

/** Allocating child node.
  *
  * Given a child index, this routine creates an empty ProjectedNode child node
  * with the appropriate translation and Hilbert path parameters. */
template<int D>
void ProjectedNode<D>::createChild(int cIdx) {
    NOT_IMPLEMENTED_ABORT;
//    assert(this->children != 0);
//    assert(this->children[cIdx] == 0);
//    ProjectedNode<D> *child = new ProjectedNode<D>(*this, cIdx);
//    this->children[cIdx] = child;
}

/** Generating child node.
  *
  * This routine creates 2^D children empty GenNodes with the appropriate
  * translation and Hilbert path parameters. */
template<int D>
void ProjectedNode<D>::genChild(int cIdx) {
    assert(this->children[cIdx] == 0);
    MRNode<D> *child = new GenNode<D>(*this, cIdx);
    this->children[cIdx] = child;
}

template class ProjectedNode<1> ;
template class ProjectedNode<2> ;
template class ProjectedNode<3> ;

