// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>

#ifndef CGAL_STRAIGHT_SKELETON_2_H
#define CGAL_STRAIGHT_SKELETON_2_H 1

#ifndef CGAL_STRAIGHT_SKELETON_ITEMS_2_H
#include <CGAL/Straight_skeleton_items_2.h>
#endif

#ifndef CGAL_HALFEDGEDS_DEFAULT_H
#include <CGAL/HalfedgeDS_default.h>
#endif

CGAL_BEGIN_NAMESPACE

template<  class Traits_
         , class Items_ = Straight_skeleton_items_2
         , class Alloc_ = CGAL_ALLOCATOR(int)
        >
class Straight_skeleton_2 : public CGAL_HALFEDGEDS_DEFAULT <Traits_,Items_,Alloc_>
{
public :

  typedef Traits_ Traits ;

  typedef Straight_skeleton_2<Traits_,Items_,Alloc_> Self ;
  
  typedef CGAL_HALFEDGEDS_DEFAULT <Traits_,Items_,Alloc_> Base ;
  
  typedef typename Base::Vertex_base     Vertex ;
  typedef typename Base::Halfedge_base   Halfedge ;
  typedef typename Base::Face_base       Face ;
  
  typedef typename Base::Vertex_handle   Vertex_handle ;
  typedef typename Base::Halfedge_handle Halfedge_handle  ;
  typedef typename Base::Face_handle     Face_handle  ;
  
  typedef typename Base::Vertex_iterator Vertex_iterator ;
  typedef typename Base::Halfedge_iterator Halfedge_iterator  ;
  typedef typename Base::Face_iterator Face_iterator  ;

  Straight_skeleton_2() {}
  
private :

    Vertex_handle vertices_push_back( const Vertex& v) { return Vertex_handle(); }
    Halfedge_handle edges_push_back( const Halfedge& h, const Halfedge& g) { return Halfedge_handle(); }
    Halfedge_handle edges_push_back( const Halfedge& h) { return Halfedge_handle(); }
    Face_handle faces_push_back( const Face& f) { return Face_handle(); }
    
    void vertices_pop_front() {}
    void vertices_pop_back() {}
    void vertices_erase( Vertex_handle v) {}
    void vertices_erase( Vertex_iterator first, Vertex_iterator last) {}
    void edges_erase( Halfedge_handle h) {}
    void edges_pop_front() {}
    void edges_pop_back()  {}
    void edges_erase( Halfedge_iterator first, Halfedge_iterator last) {}
    void faces_pop_front() {}
    void faces_pop_back() {}
    void faces_erase( Face_handle f) {}
    void faces_erase( Face_iterator first, Face_iterator last) {}
    void vertices_clear() {}
    void edges_clear() {}
    void faces_clear() {}
    void clear() {}
    void vertices_splice( Vertex_iterator target, Self &source,
                          Vertex_iterator begin, Vertex_iterator end) {}

    void halfedges_splice( Halfedge_iterator target, Self &source,
                           Halfedge_iterator begin, Halfedge_iterator end) {}

    void faces_splice( Face_iterator target, Self &source,
                       Face_iterator begin, Face_iterator end) {}
    void normalize_border() {}
};


CGAL_END_NAMESPACE


#endif // CGAL_STRAIGHT_SKELETON_2_H //
// EOF //

