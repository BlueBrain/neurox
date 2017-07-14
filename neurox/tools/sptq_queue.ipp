/*
 Copyright (c) 2016, Blue Brain Project
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification,
 are permitted provided that the following conditions are met:
 1. Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software
 without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef tqueue_ipp_
#define tqueue_ipp_

namespace neurox{
namespace tools {

    // splay tree + bin queue limited to fixed step method
    // for event-sets or priority queues
    // this starts from the sptqueue.cpp file and adds a bin queue

    /* Derived from David Brower's c translation of pascal code by
     Douglas Jones.
     */
    /* The original c code is included from this file but note that instead
     of struct _spblk, we are really using TQItem
     */
    /*
     *  The following code implements the basic operations on
     *  an event-set or priority-queue implemented using splay trees:
     *
     Hines changed to void spinit(SPTREE**) for use with TQueue.
     *  SPTREE *spinit( compare )	Make a new tree
     *  SPBLK *spenq( n, q )	Insert n in q after all equal keys.
     *  SPBLK *spdeq( np )		Return first key under *np, removing it.
     *  void splay( n, q )		n (already in q) becomes the root.
     *  int n = sphead( q )         n is the head item in q (not removed).
     *
     *  In the above, n points to an SPBLK type, while q points to an
     *  SPTREE.
     *
     *  The implementation used here is based on the implementation
     *  which was used in the tests of splay trees reported in:
     *
     *    An Empirical Comparison of Priority-Queue and Event-Set Implementations,
     *	by Douglas W. Jones, Comm. ACM 29, 4 (Apr. 1986) 300-311.
     *
     *  The changes made include the addition of the enqprior
     *  operation and the addition of up-links to allow for the splay
     *  operation.  The basic splay tree algorithms were originally
     *  presented in:
     *
     *	Self Adjusting Binary Trees,
     *		by D. D. Sleator and R. E. Tarjan,
     *			Proc. ACM SIGACT Symposium on Theory
     *			of Computing (Boston, Apr 1983) 235-245.
     *
     *  The enq and enqprior routines use variations on the
     *  top-down splay operation, while the splay routine is bottom-up.
     *  All are coded for speed.
     *
     *  Written by:
     *    Douglas W. Jones
     *
     *  Translated to C by:
     *    David Brower, daveb@rtech.uucp
     *
     * Thu Oct  6 12:11:33 PDT 1988 (daveb) Fixed spdeq, which was broken
     *	handling one-sptq_node trees.  I botched the pascal translation of
     *	a VAR parameter.
     *
     * July 2016 - Tim Ewart - timothee.ewart@gmail.com
     *  - namespace stdq instead of C style signature
     *  - Remove all macros
     *  - Item becomes sptq_node
     *  - Introduce generacity for the value of the sptq_node and integrate the
     template argumenent comparator
     *    to respect the std::API
     *  - Introduce object function to select the good comparator for std::less and
     std::greater
     * Sep. 2016 - Tim Ewart - timothee.ewart@gmail.com
     *  - simplify compartor function
     *  - replug spdelete for move function
     */

    /** comparator */
    template <class T, class Compare>
    inline bool comparator(const T& a, const T& b, Compare c = Compare()) {
        return c(a, b);
    }

    /*----------------
     *
     * spinit() -- initialize an empty splay tree
     *
     */
    template <class T>
    void spinit(SPTREE<T>* q) {
        q->root = NULL;
    }

    /*----------------
     *
     *  spenq() -- insert item in a tree.
     *
     *  put n in q after all other nodes with the same key; when this is
     *  done, n will be the root of the splay tree representing q, all nodes
     *  in q with keys less than or equal to that of n will be in the
     *  left subtree, all with greater keys will be in the right subtree;
     *  the tree is split into these subtrees from the top down, with rotations
     *  performed along the way to shorten the left branch of the right subtree
     *  and the right branch of the left subtree
     */
    template <class T, class Compare>
    sptq_node<T>* spenq(sptq_node<T>* n, SPTREE<T>* q) {
        typename sptq_node<T>::value_type key;

        sptq_node<T>* left;  /* the rightmost sptq_node in the left tree */
        sptq_node<T>* right; /* the leftmost sptq_node in the right tree */
        sptq_node<T>* next;  /* the root of the unsplit part */
        sptq_node<T>* temp;

        n->parent_ = NULL;
        next = q->root;
        q->root = n;
        if (next == NULL) /* trivial enq */
        {
            n->left_ = NULL;
            n->right_ = NULL;
        } else /* difficult enq */
        {
            key = n->t_;
            left = n;
            right = n;

            /* n's left and right children will hold the right and left
             splayed trees resulting from splitting on n->key;
             note that the children will be reversed! */

            if (comparator<T, Compare>(next->t_, key))
                goto two;

        one: /* assert next->key <= key */

            do /* walk to the right in the left tree */
            {
                temp = next->right_;
                if (temp == NULL) {
                    left->right_ = next;
                    next->parent_ = left;
                    right->left_ = NULL;
                    goto done; /* job done, entire tree split */
                }

                if (comparator<T, Compare>(temp->t_, key)) {
                    left->right_ = next;
                    next->parent_ = left;
                    left = next;
                    next = temp;
                    goto two; /* change sides */
                }

                next->right_ = temp->left_;
                if (temp->left_ != NULL)
                    temp->left_->parent_ = next;
                left->right_ = temp;
                temp->parent_ = left;
                temp->left_ = next;
                next->parent_ = temp;
                left = temp;
                next = temp->right_;
                if (next == NULL) {
                    right->left_ = NULL;
                    goto done; /* job done, entire tree split */
                }

            } while (!comparator<T, Compare>(next->t_, key)); /* change sides */

        two: /* assert next->key > key */

            do /* walk to the left in the right tree */
            {
                temp = next->left_;
                if (temp == NULL) {
                    right->left_ = next;
                    next->parent_ = right;
                    left->right_ = NULL;
                    goto done; /* job done, entire tree split */
                }

                if (!comparator<T, Compare>(temp->t_, key)) {
                    right->left_ = next;
                    next->parent_ = right;
                    right = next;
                    next = temp;
                    goto one; /* change sides */
                }
                next->left_ = temp->right_;
                if (temp->right_ != NULL)
                    temp->right_->parent_ = next;
                right->left_ = temp;
                temp->parent_ = right;
                temp->right_ = next;
                next->parent_ = temp;
                right = temp;
                next = temp->left_;
                if (next == NULL) {
                    left->right_ = NULL;
                    goto done; /* job done, entire tree split */
                }

            } while (comparator<T, Compare>(next->t_, key)); /* change sides */

            goto one;

        done: /* split is done, branches of n need reversal */

            temp = n->left_;
            n->left_ = n->right_;
            n->right_ = temp;
        }

        return (n);

    } /* spenq */

    /*----------------
     *
     *  spdeq() -- return and remove head sptq_node from a subtree.
     *
     *  remove and return the head sptq_node from the sptq_node set; this deletes
     *  (and returns) the leftmost sptq_node from q, replacing it with its right
     *  subtree (if there is one); on the way to the leftmost sptq_node, rotations
     *  are performed to shorten the left branch of the tree
     */
    template <class T>
    sptq_node<T>* spdeq(sptq_node<T>** np) /* pointer to a sptq_node pointer */
    {
        sptq_node<T>* deq;        /* one to return */
        sptq_node<T>* next;       /* the next thing to deal with */
        sptq_node<T>* left;       /* the left child of next */
        sptq_node<T>* farleft;    /* the left child of left */
        sptq_node<T>* farfarleft; /* the left child of farleft */

        if (np == NULL || *np == NULL) {
            deq = NULL;
        } else {
            next = *np;
            left = next->left_;
            if (left == NULL) {
                deq = next;
                *np = next->right_;

                if (*np != NULL)
                    (*np)->parent_ = NULL;

            } else
                for (;;) /* left is not null */
                {
                    /* next is not it, left is not NULL, might be it */
                    farleft = left->left_;
                    if (farleft == NULL) {
                        deq = left;
                        next->left_ = left->right_;
                        if (left->right_ != NULL)
                            left->right_->parent_ = next;
                        break;
                    }

                    /* next, left are not it, farleft is not NULL, might be it */
                    farfarleft = farleft->left_;
                    if (farfarleft == NULL) {
                        deq = farleft;
                        left->left_ = farleft->right_;
                        if (farleft->right_ != NULL)
                            farleft->right_->parent_ = left;
                        break;
                    }

                    /* next, left, farleft are not it, rotate */
                    next->left_ = farleft;
                    farleft->parent_ = next;
                    left->left_ = farleft->right_;
                    if (farleft->right_ != NULL)
                        farleft->right_->parent_ = left;
                    farleft->right_ = left;
                    left->parent_ = farleft;
                    next = farleft;
                    left = farfarleft;
                }
        }

        return (deq);

    } /* spdeq */

    /*----------------
     *
     *  splay() -- reorganize the tree.
     *
     *  the tree is reorganized so that n is the root of the
     *  splay tree representing q; results are unpredictable if n is not
     *  in q to start with; q is split from n up to the old root, with all
     *  nodes to the left of n ending up in the left subtree, and all nodes
     *  to the right of n ending up in the right subtree; the left branch of
     *  the right subtree and the right branch of the left subtree are
     *  shortened in the process
     *
     *  this code assumes that n is not NULL and is in q; it can sometimes
     *  detect n not in q and complain
     */
    template <class T>
    void splay(sptq_node<T>* n, SPTREE<T>* q) {
        sptq_node<T>* up;     /* points to the sptq_node being dealt with */
        sptq_node<T>* prev;   /* a descendent of up, already dealt with */
        sptq_node<T>* upup;   /* the parent of up */
        sptq_node<T>* upupup; /* the grandparent of up */
        sptq_node<T>* left;   /* the top of left subtree being built */
        sptq_node<T>* right;  /* the top of right subtree being built */

        left = n->left_;
        right = n->right_;
        prev = n;
        up = prev->parent_;

        while (up != NULL) {
            /* walk up the tree towards the root, splaying all to the left of
             n into the left subtree, all to right into the right subtree */

            upup = up->parent_;
            if (up->left_ == prev) /* up is to the right of n */
            {
                if (upup != NULL && upup->left_ == up) /* rotate */
                {
                    upupup = upup->parent_;
                    upup->left_ = up->right_;
                    if (upup->left_ != NULL)
                        upup->left_->parent_ = upup;
                    up->right_ = upup;
                    upup->parent_ = up;
                    if (upupup == NULL)
                        q->root = up;
                    else if (upupup->left_ == upup)
                        upupup->left_ = up;
                    else
                        upupup->right_ = up;
                    up->parent_ = upupup;
                    upup = upupup;
                }
                up->left_ = right;
                if (right != NULL)
                    right->parent_ = up;
                right = up;

            } else /* up is to the left of n */
            {
                if (upup != NULL && upup->right_ == up) /* rotate */
                { /** (Tim's comment), Please Mario leaves the game, up up up big up is
                     not a new life */
                    upupup = upup->parent_;
                    upup->right_ = up->left_;
                    if (upup->right_ != NULL)
                        upup->right_->parent_ = upup;
                    up->left_ = upup;
                    upup->parent_ = up;
                    if (upupup == NULL)
                        q->root = up;
                    else if (upupup->right_ == upup)
                        upupup->right_ = up;
                    else
                        upupup->left_ = up;
                    up->parent_ = upupup;
                    upup = upupup;
                }
                up->right_ = left;
                if (left != NULL)
                    left->parent_ = up;
                left = up;
            }
            prev = up;
            up = upup;
        }

#ifdef DEBUG
        if (q->root != prev) {
            /*	fprintf(stderr, " *** bug in splay: n not in q *** " ); */
            abort();
        }
#endif

        n->left_ = left;
        n->right_ = right;
        if (left != NULL)
            left->parent_ = n;
        if (right != NULL)
            right->parent_ = n;
        q->root = n;
        n->parent_ = NULL;

    } /* splay */

    /*----------------
     *
     * sphead() --  return the "lowest" element in the tree.
     *
     *      returns a reference to the head event in the event-set q,
     *      represented as a splay tree; q->root ends up pointing to the head
     *      event, and the old left branch of q is shortened, as if q had
     *      been splayed about the head element; this is done by dequeueing
     *      the head and then making the resulting queue the right son of
     *      the head returned by spdeq; an alternative is provided which
     *      avoids splaying but just searches for and returns a pointer to
     *      the bottom of the left branch
     */
    template <class T>
    sptq_node<T>* sphead(SPTREE<T>* q) {
        sptq_node<T>* x;

        /* splay version, good amortized bound */
        x = spdeq(&q->root);
        if (x != NULL) {
            x->right_ = q->root;
            x->left_ = NULL;
            x->parent_ = NULL;
            if (q->root != NULL)
                q->root->parent_ = x;
        }
        q->root = x;

        /* alternative version, bad amortized bound,
         but faster on the average */

        return (x);

    } /* sphead */

    /*----------------
     *
     * spdelete() -- Delete sptq_node from a tree.
     *
     *	n is deleted from q; the resulting splay tree has been splayed
     *	around its new root, which is the successor of n
     *
     */
    template <class T>
    void spdelete(sptq_node<T>* n, SPTREE<T>* q) {
        sptq_node<T>* x;

        splay(n, q);
        x = spdeq(&q->root->right_);
        if (x == NULL) /* empty right subtree */
        {
            q->root = q->root->left_;
            if (q->root)
                q->root->parent_ = NULL;
        } else /* non-empty right subtree */
        {
            x->parent_ = NULL;
            x->left_ = q->root->left_;
            x->right_ = q->root->right_;
            if (x->left_ != NULL)
                x->left_->parent_ = x;
            if (x->right_ != NULL)
                x->right_->parent_ = x;
            q->root = x;
        }
    } /* spdelete */

} //tools
} //neurox

#endif
