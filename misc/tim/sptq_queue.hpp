/*
Copyright (c) 2016, Blue Brain Project
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
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

#ifndef tqueue_hpp_
#define tqueue_hpp_

#include <assert.h>
#include <cstring>
#include <functional>
#include <ostream>
#include <stdio.h>

#include "node.h"

namespace tim {
    /** The queue: TQeue from Michael starts here, not compliant with std for the
    container,
    but ok for the type support and the comparator, by default std::less and not
    std::greated
    as it was in the original version, already a bit of clean up specially for the
    push ... */

    // node for the splay tree
    template <class T>
    struct sptq_node : TQitem {  // node is public because struct
        typedef T value_type;
        explicit sptq_node(value_type t = value_type()) : t_(t), left_(0), right_(0), parent_(0){};
        value_type t_; /*!< Data of the node */
        sptq_node* left_;
        sptq_node* right_;
        sptq_node* parent_;
    };

    template <class T>
    struct SPTREE {
        sptq_node<T>* root; /* root node */
    };

    /** Forward declarations for the c++ interface */
    /* init tree */
    template <class T>
    void spinit(SPTREE<T>*);

    /* insert item into the tree */
    template <class T, class Compare>
    sptq_node<T>* spenq(sptq_node<T>*, SPTREE<T>*);

    /* return and remove lowest item in subtree */
    template <class T>
    sptq_node<T>* spdeq(sptq_node<T>**);

    /* return first node in tree */
    template <class T>
    sptq_node<T>* sphead(SPTREE<T>*);

    /*return a looking node */
    template <class T>
    void spdelete(sptq_node<T>*, SPTREE<T>*);

    template <class T, class Compare = std::less<T> >
    class sptq_queue {
      public:
        typedef SPTREE<T> container;
        typedef std::size_t size_type;
        typedef T value_type;
        typedef sptq_node<T> node_type;

        /** \fn sptq_queue()
         *  \brief initializes the sptq queue
         */
        inline sptq_queue() : size_(0) {
            spinit(&q);
        }

        /** \fn ~sptq_queue()
         *  \brief destructor clean up the queue using clear() function
         */
        inline ~sptq_queue() {
            clear();
        }

        /** \fn clear()
         *  \brief destructor delete all node
         */
        inline void clear() {
            node_type* n;
            while ((n = spdeq(&(&q)->root)) != NULL)
                delete n;
        }

        /** \fn push(const value_type& t)
         *  \brief push a element into the queue, allocate a node
         *  \param t the value of the element to push
         */
        inline void push(const value_type& value) {
            node_type* n = new node_type(value);
            spenq<T, Compare>(n, &q);  // the Comparator is use only here
            size_++;
        }

        /** \fn push(node_type* t)
         *  \brief push a node into the queue, the node must be allocated BEFORE
         *  \param n the node to push
         */
        inline void push(node_type* n) {
            spenq<T, Compare>(n, &q);
            size_++;
        }

        /** \fn pop()
         *  \brief Removes the top element from the priority queue
         */
        inline void pop() {
            if (!empty()) {
                node_type* n = spdeq(&(&q)->root);
                delete n;  // pop remove definitively the element else memory leak
                size_--;
            }
        }

        /** \fn top()
         *  \brief Return the top element from the priority queue
         *  \return Return the top element from the priority queue
         */
        inline value_type top() {
            value_type tmp = value_type();
            if (!empty())
                tmp = sphead<T>(&q)->t_;
            return tmp;
        }

        /** \fn size()
         *  \brief returns the number of elements
         */
        inline size_type size() {
            return size_;
        }

        /** \fn empty()
         *  \brief checks whether the underlying container is empty
         *  \return true if empty else false
         */
        inline bool empty() {
            return !bool(size_);  // is it true on Power?
        }

        /** \fn find(node_type* n)
         *  \brief Find and remove a looking node into the priority_queue
         *  \param n The looking node
         */
        inline node_type* find(node_type* n) {
            spdelete(n, &q);
            size_--;  // WARNING remove the node but do not delete it
            return n;
        }

      private:
        size_type size_; /*!< size of the queue */
        container q;     /*!< the legacy "C" data structure to store the node */
    };
}  // end namespace
#include "sptq_queue.ipp"

#endif
