// Copyright (c) 2024 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef HEXAHEDRAL_SUBDIVISION_H
#define HEXAHEDRAL_SUBDIVISION_H

#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_basic_operations.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Unique_hash_map.h>

#include "Logger.h"
#include "Volume_info_for_hexahedral_subdivision.h"

using namespace std;

typedef unsigned int level_size;

/// Return log2(x) if x is a power of 2.
level_size integer_log2(level_size x)
{ return __builtin_ctz(x); }

/// Test if x is a power of 2 (i.e. 2^1, ..., 2^32).
///   (note that this function return false for x==0)
bool is_power_2(level_size x)
{ return (__builtin_popcount(x)==1); }

/// @return true iff x is even
bool is_even(level_size x)
{
  // return (__builtin_parity(x)==0);
  return (x & 1)==0;
}

// Class to use fractions, denominator always being a power of 2.
// Fraction are always simplified (i.e. there is no 2 multiplicator in denominator)
// @TODO change the simplification and assertions if we want to define other subdvision scheme
//       where we do not always cut in 2.
class CFraction
{
public:
  // Contruct a fraction; den must be a power of 2.
  CFraction(level_size num, level_size den): m_num(num), m_den(den)
  {
    assert(is_power_2(den));
    simplify();
  }

  // Add two fractions (this and 1/e); result is stored in this.
  CFraction& operator+=(level_size e)
  {
    if(m_den!=e)
    {
      m_num = m_num * e + m_den;
      m_den *= e;
    }
    else
    {
      ++m_num;
    }
    simplify();

    assert(is_power_2(m_den));

    return *this;
  }

  CFraction& operator/=(level_size e)
  {
    m_den *= e;
    return *this;
  }

  // Add two fractions (this and f); result is stored in this.
  CFraction& operator+=(const CFraction& f)
  {
    if(m_den!=f.m_den)
    {
      m_num = m_num * f.m_den + f.m_num * m_den;
      m_den *= f.m_den;
    }
    else
    {
      m_num += f.m_num;
    }
    simplify();

    assert(is_power_2(m_den));

    return *this;
  }

  bool operator==(const CFraction& f) const
  {
    // Here this works because fractions are always simplified.
    return m_num==f.m_num && m_den==f.m_den;
  }

  bool operator!=(const CFraction& f) const
  { return !operator==(f); }

  bool operator<(const CFraction& f) const
  {
    if (m_den==f.m_den) return m_num<f.m_num;
    return (m_num*f.m_den)<(f.m_num*m_den);
  }

  bool operator<=(const CFraction& f) const
  {
    if (m_den==f.m_den) return m_num<=f.m_num;
    return (m_num*f.m_den)<=(f.m_num*m_den);
  }

  bool operator>(const CFraction& f) const
  { return f.operator<(*this); }

  bool operator>=(const CFraction& f) const
  { return f.operator<=(*this); }

  level_size get_num() const
  { return m_num; }

  level_size get_den() const
  { return m_den; }

protected:
  /// Simplify the fraction.
  void simplify()
  {
    assert(is_power_2(m_den));
    while(is_even(m_num) && is_even(m_den))
    {
      m_num>>=1;
      m_den>>=1;
    }
  }

  level_size m_num, m_den;
};

template<typename LCC>
class Hexahedral_subdivision
{
public:
  typedef Hexahedral_subdivision<LCC> Self;
  typedef typename LCC::Dart_descriptor Dart_descriptor;
  typedef typename LCC::Vertex_attribute_handle Vertex_attribute_handle;
  typedef typename LCC::Point Point;
  typedef Volume_info_for_hexahedral_subdivision<LCC> Volume_info;

  // Default constructor: the Linear_cell_complex_subdivision is not initialized.
  Hexahedral_subdivision() : m_lcc(nullptr),
    m_split_edge_in_two(nullptr),
    m_split_edge_in_two_old(nullptr),
    m_split_face_in_four(nullptr),
    m_split_face_in_four_old(nullptr)
  {}

  // Constructor taking a lcc as parameter: the Linear_cell_complex_subdivision is initialized.
  Hexahedral_subdivision(LCC * alcc)  : m_lcc(nullptr),
    m_split_edge_in_two(nullptr),
    m_split_edge_in_two_old(nullptr),
    m_split_face_in_four(nullptr),
    m_split_face_in_four_old(nullptr)
  { init(alcc); }

  Volume_info& volume_info(typename LCC::template Attribute_handle<3>::type ah)
  {
    // TODO template specialisation 
#if 1 // WITH_INDEX
    return m_volume_info[ah];
#else
    return m_volume_info[ah->id()];
#endif
  }
  const Volume_info& volume_info(typename LCC::template Attribute_handle<3>::type ah) const
  { return m_volume_info[ah->id()]; }

  Volume_info& volume_info(Dart_descriptor dh)
  { return volume_info(m_lcc->template attribute<3>(dh)); }
  const Volume_info& volume_info(Dart_descriptor dh) const
  { return volume_info(m_lcc->template attribute<3>(dh)); }

  // Init: we reserve a mark for corner darts; and initialize all
  // the subdivision levels to 1.
  // Moreover we create a Volume_lcc and associate it to each hexahedron
  // (this allows to process mixed meshes while allowing to subdivide hexahedra).
  // @pre takes in parameter a linear cell complex which is not yet subdivided.
  void init(LCC * alcc)
  {
    assert(m_lcc==nullptr);
    m_lcc = alcc;
    m_corner_mark = m_lcc->get_new_mark();
    m_edge_first_direction = m_lcc->get_new_mark();

    m_volume_info.resize(m_lcc->template attributes<3>().capacity());

    typename LCC::size_type test_mark = m_lcc->get_new_mark();

    for(typename LCC::Dart_range::iterator it(m_lcc->darts().begin()),
          itend(m_lcc->darts().end()) ; it!=itend ; ++it)
    {
      if(!m_lcc->is_marked(it, test_mark))
      {
        if(m_lcc->is_volume_combinatorial_hexahedron(it))
        {
          Dart_descriptor mindh=it;
          assert(m_lcc->template attribute<3>(it)!=m_lcc->null_handle);

          for (typename LCC::template Dart_of_cell_basic_range<3>::iterator
                 itvol   =m_lcc->template darts_of_cell_basic<3>(it, test_mark).begin(),
                 itvolend=m_lcc->template darts_of_cell_basic<3>(it, test_mark).end();
               itvol!=itvolend ; ++itvol)
          {
            m_lcc->mark(itvol, test_mark);
            m_lcc->mark(itvol, m_corner_mark);
            // Not necessary: 3-attributes are already created
            //        m_lcc->template set_dart_attribute<3>(itvol, a);
            if ( m_lcc->point(itvol)<m_lcc->point(mindh) ||
                 (m_lcc->point(itvol)==m_lcc->point(mindh) &&
                  m_lcc->point(m_lcc->other_extremity(itvol))<
                  m_lcc->point(m_lcc->other_extremity(mindh))) )
              mindh=itvol;
          }
          volume_info(mindh).initialize_corners(*m_lcc, mindh);
          mark_hexahedron(mindh);
        }
        else
        {
          CGAL::mark_cell<LCC,3>(*m_lcc, it, test_mark);
        }
      }
    }

    m_lcc->free_mark(test_mark);
  }

  void init_one_hexa(Dart_descriptor it)
  {
    /* typename LCC::template Attribute_handle<3>::type
        a=m_lcc->template create_attribute<3>(); */

    Dart_descriptor mindh=it;
    for (typename LCC::template Dart_of_cell_range<3>::iterator
           itvol   =m_lcc->template darts_of_cell<3>(it).begin(),
           itvolend=m_lcc->template darts_of_cell<3>(it).end();
         itvol!=itvolend ; ++itvol)
    {
      m_lcc->mark(itvol, m_corner_mark);
      // m_lcc->template set_dart_attribute<3>(itvol, a);
      if ( m_lcc->point(itvol)<m_lcc->point(mindh) ||
           (m_lcc->point(itvol)==m_lcc->point(mindh) &&
            m_lcc->point(m_lcc->other_extremity(itvol))<
            m_lcc->point(m_lcc->other_extremity(mindh))) )
        mindh=itvol;
    }
    volume_info(mindh).initialize_corners(*m_lcc, mindh);
    mark_hexahedron(mindh);
  }

  void init_without_3attribs(LCC * alcc)
  {
    assert(m_lcc==nullptr);
    m_lcc = alcc;
    m_corner_mark = m_lcc->get_new_mark();
    m_edge_first_direction = m_lcc->get_new_mark();

    typename LCC::size_type test_mark = m_lcc->get_new_mark();

    for(typename LCC::Dart_range::iterator it(m_lcc->darts().begin()),
          itend(m_lcc->darts().end()) ; it!=itend ; ++it)
    {
      if(!m_lcc->is_marked(it, test_mark))
      {
        if(m_lcc->is_face_combinatorial_polygon(it, 4)) // A square
        {
          set_first_direction(it);
          set_first_direction(m_lcc->template beta<1, 1>(it));
          m_lcc->mark(it, m_corner_mark);
          m_lcc->mark(m_lcc->template beta<1>(it), m_corner_mark);
          m_lcc->mark(m_lcc->template beta<1, 1>(it), m_corner_mark);
          m_lcc->mark(m_lcc->template beta<0>(it), m_corner_mark);

          if ( !m_lcc->template is_free<3>(it) )
          {
            set_first_direction(m_lcc->template beta<3>(it));
            set_first_direction(m_lcc->template beta<1, 1, 3>(it));
            m_lcc->mark(m_lcc->template beta<3>(it), m_corner_mark);
            m_lcc->mark(m_lcc->template beta<1, 3>(it), m_corner_mark);
            m_lcc->mark(m_lcc->template beta<1, 1, 3>(it), m_corner_mark);
            m_lcc->mark(m_lcc->template beta<0, 3>(it), m_corner_mark);
          }
        }
        CGAL::mark_cell<LCC, 2>(*m_lcc, it, test_mark);
      }
    }
    m_lcc->free_mark(test_mark);
  }

  // Destructor
  ~Hexahedral_subdivision()
  {
    if (m_lcc!=nullptr)
    {
      m_lcc->free_mark(m_corner_mark);
      m_lcc->free_mark(m_edge_first_direction);
    }
  }

  LCC* get_lcc()
  { return m_lcc; }

  // Debug function
  void debug(const char* txt="")
  {
    if (txt[0]!=0)
    {
      std::cout<<"****************************"<<txt<<"****************************"<<std::endl;
    }

    std::cout<<"Number of hexahedra: "<<m_lcc->template number_of_attributes<3>()<<std::endl;

    std::cout<<"Number of darts: "<<m_lcc->number_of_darts()<<std::endl
            <<"Number of corners: "<<m_lcc->number_of_marked_darts(m_corner_mark)<<std::endl
           <<"Number of darts first dir: "<<m_lcc->number_of_marked_darts(m_edge_first_direction)<<std::endl
          <<"Number of darts second dir: "<<m_lcc->number_of_darts()-m_lcc->number_of_marked_darts(m_edge_first_direction)<<std::endl;

    m_lcc->display_characteristics(std::cout);
    std::cout<<", valid="<<m_lcc->is_valid()<<std::endl;
  }

  bool check_corner_darts()
  {
    typename LCC::size_type m=m_lcc->get_new_mark();

    for (typename LCC::template Attribute_range<3>::type::iterator
         it=m_lcc->template attributes<3>().begin(),
         itend=m_lcc->template attributes<3>().end();
         it!=itend; ++it)
    {
      // First test that darts in corner array are marked
      for (int i=0; i<6; ++i)
        for (int j=0; j<4; ++j)
        {
          if (!is_corner_dart(volume_info(it).get_corner_dart(i, j)))
          {
            std::cout<<"Problem: a corner dart is not marked corner"<<std::endl;
            m_lcc->free_mark(m);
            return false;
          }
          else
          {
            m_lcc->mark(volume_info(it).get_corner_dart(i, j), m);
          }
          if (!is_hierarchical_edge(volume_info(it).get_corner_dart(i, j),
                                    volume_info(it).get_corner_dart(i, (j+1)%4)))
          {
            std::cout<<"Problem: edge "<<m_lcc->point(volume_info(it).get_corner_dart(i, j))
                    <<" -> "<<m_lcc->point(volume_info(it).get_corner_dart(i, (j+1)%4))
                   <<" is not correct."<<std::endl;
            display_hierarchical_edge(volume_info(it).get_corner_dart(i, j),
                                      volume_info(it).get_corner_dart(i, (j+1)%4));
            m_lcc->free_mark(m);
            return false;
          }
        }
      // Second test that other darts are not marked
      for (typename LCC::template Dart_of_cell_range<3>::iterator
           it2=m_lcc->template darts_of_cell<3>(it->dart()).begin(),
           it2end=m_lcc->template darts_of_cell<3>(it->dart()).end();
           it2!=it2end; ++it2)
      {
        if ( !m_lcc->is_marked(it2, m) && is_corner_dart(it2))
        {
          std::cout<<"Problem: a non corner dart is marked corner"<<std::endl;
          m_lcc->free_mark(m);
          return false;
        }
        m_lcc->unmark(it2, m);
      }
    }

    m_lcc->free_mark(m);
    return true;
  }

  void display_corner_darts(Dart_descriptor (&cd) [6][4])
  {
    std::cout<<"********* Corners darts *********"<<std::endl;
    for (int i=0; i<6; ++i)
    {
      std::cout<<"Face "<<i<<": ";
      for (int j=0; j<4; ++j)
      {
        std::cout<<m_lcc->point(cd[i][j])<<"  ";
      }
      std::cout<<std::endl;
    }
  }
  void display_middle_darts(Dart_descriptor (&md) [6][4])
  {
    std::cout<<"********* Middle darts *********"<<std::endl;
    for (int i=0; i<6; ++i)
    {
      std::cout<<"Face "<<i<<": ";
      for (int j=0; j<4; ++j)
      {
        std::cout<<m_lcc->point(md[i][j])<<"  ";
      }
      std::cout<<std::endl;
    }
  }
  void display_face_darts(Dart_descriptor (&fd) [6])
  {
    std::cout<<"********* Face darts *********"<<std::endl;
    for (int i=0; i<6; ++i)
    {
      std::cout<<"Face "<<i<<": ";
      std::cout<<m_lcc->point(fd[i])<<std::endl;
    }
  }

  /** Compute the length of the border of a face (a border is the sequence
   *  of darts of the face with same direction).
   *  @param out the subdvision level of the edges of the same border
   *  @return the first dart of the next border
   */
  Dart_descriptor compute_length_of_face_border(Dart_descriptor dh, level_size& level,
                                            level_size& oppositelevel)
  {
    assert(is_first_direction(dh)!=
        is_first_direction(m_lcc->template beta<0>(dh)));

    CFraction sum(0, 1);
    CFraction sum2(0, 1);
    Dart_descriptor cur=dh;

    while ( is_first_direction(cur)==is_first_direction(dh) )
    {
      sum += get_subdivision_level_of_edge(cur);
      if ( !m_lcc->template is_free<3>(cur))
        sum2 += get_subdivision_level_of_edge(m_lcc->template beta<3>(cur));

      cur=m_lcc->template beta<1>(cur);
    }
    assert(sum.get_num()==1);
    level=sum.get_den();
    oppositelevel=sum2.get_den();

    return cur;
  }

  bool is_hierarchical_edge(Dart_descriptor dh1, Dart_descriptor dh2)
  {
    CFraction size_he(1, 1);
    CFraction sum(0, 1);
    while(sum<size_he)
    {
      sum += get_subdivision_level_of_edge(dh1);
      if (sum<size_he) dh1 = get_next_edge_dart(dh1);
      else dh1 = m_lcc->template beta<1>(dh1);
    }
    assert(sum==size_he);
    return (dh1==dh2);
  }

  void display_hexahedron(Dart_descriptor dh)
  {
    typename LCC::size_type m=m_lcc->get_new_mark();

    for(typename LCC::template Dart_of_cell_range<3>::iterator
          it(m_lcc->template darts_of_cell<3>(dh).begin()),
          itend(m_lcc->template darts_of_cell<3>(dh).end());
        it!=itend; ++it)
    {
      if (!m_lcc->is_marked(it, m))
      {
        std::cout<<"*** Face *** "<<std::endl;
        for(typename LCC::template Dart_of_orbit_range<1>::iterator
              it2(m_lcc->template darts_of_orbit<1>(it).begin()),
              it2end(m_lcc->template darts_of_orbit<1>(it).end());
            it2!=it2end; ++it2)
        {
          std::cout<<m_lcc->point(it2)<<" dir="<<is_first_direction(it2)
                  <<" corner="<<is_corner_dart(it2)
                 <<" subdvision="<<get_subdivision_level_of_edge(it2)<<std::endl;
          m_lcc->mark(it2, m);
        }
      }
    }

    m_lcc->free_mark(m);
  }

  void display_hierarchical_edge(Dart_descriptor dh1, Dart_descriptor dh2)
  {
    CFraction size_he(1, 1);
    CFraction sum(0, 1);
    while(sum<size_he)
    {
      std::cout<<m_lcc->point(dh1)<<" dir="<<is_first_direction(dh1)
              <<" subdvidision="<<get_subdivision_level_of_edge(dh1)<<std::endl;
      sum += get_subdivision_level_of_edge(dh1);
      if (sum<size_he) dh1 = get_next_edge_dart(dh1);
      else dh1 = m_lcc->template beta<1>(dh1);
    }
    std::cout<<m_lcc->point(dh1)<<" dir="<<is_first_direction(dh1)
            <<" subdvidision="<<get_subdivision_level_of_edge(dh1)<<std::endl;
    assert(sum==size_he);
    if (dh2!=m_lcc->null_handle)
      std::cout<<"Next corner="<<m_lcc->point(dh2)<<std::endl;
  }

  void display_hierarchical_edge(Dart_descriptor dh1)
  { display_hierarchical_edge(dh1, m_lcc->null_handle); }

  // Subdivide the 12 edges of a given hexahedron which are not yet subdivided.
  void subdivide_edges_of_hexahedron(Dart_descriptor dh)
  {
    // Logger logger("Linear_cell_complex_subdivision::subdivide_edges_of_hexahedron");

    assert(m_lcc->template attribute<3>(dh)!=m_lcc->null_handle);

    for(int i=0; i<12; i++)
    {
      if (!is_edge_already_subdivided(volume_info(dh).get_edge(i)))
        split_edge_in_two(volume_info(dh).get_edge(i));
    }
  }

  // Split the ith face of a given hexahedron in 4.
  void split_one_face_of_hexahedron_in_four(Dart_descriptor dh, unsigned int face)
  {
    assert(m_lcc->template attribute<3>(dh)!=m_lcc->null_handle);

    for (int i=0; i<4; ++i)
    {
      if (!is_edge_already_subdivided(m_lcc->template attribute<3>(dh)->info().
                                      get_corner_dart(face, i)))
        split_edge_in_two(m_lcc->template attribute<3>(dh)->info().
                          get_corner_dart(face, i));
    }

    if (!is_face_already_subdivided(m_lcc->template attribute<3>(dh)->info().
                                    get_corner_dart(face, 0)))
    {
      split_face_in_four(m_lcc->template attribute<3>(dh)->info().
                         get_corner_dart(face, 0),
                         m_lcc->template attribute<3>(dh)->info().
                         get_corner_dart(face, 1),
                         m_lcc->template attribute<3>(dh)->info().
                         get_corner_dart(face, 2),
                         m_lcc->template attribute<3>(dh)->info().
                         get_corner_dart(face, 3));
    }
  }

  // Subdivide the 6 faces of a given hexahedron which are not yet subdivided.
  // @pre all the 12 edges of the hexahedron must were already subdivided.
  void subdivide_faces_of_hexahedron(Dart_descriptor dh)
  {
    // Logger logger("Linear_cell_complex_subdivision::subdivide_faces_of_hexahedron");

    assert(m_lcc->template attribute<3>(dh)!=m_lcc->null_handle);

    for(int i=0; i<6; i++)
    {
      if (!is_face_already_subdivided(volume_info(dh).get_corner_dart(i, 0)))
        split_face_in_four(volume_info(dh).get_corner_dart(i, 0),
                           volume_info(dh).get_corner_dart(i, 1),
                           volume_info(dh).get_corner_dart(i, 2),
                           volume_info(dh).get_corner_dart(i, 3));
    }
  }

  // Subdivide a given hexahedron.
  // @pre all the 12 edges and the 12 faces of the hexahedron must were already subdivided.
  void subdivide_hexahedron(Dart_descriptor dh)
  {
    // Logger logger("Linear_cell_complex_subdivision::subdivide_hexahedron");

    assert(m_lcc->template attribute<3>(dh)!=m_lcc->null_handle);

    Dart_descriptor middledarts[6][4];
    Dart_descriptor facedarts[6];

    for(int i=0; i<6; i++)
    {
      assert(is_face_already_subdivided(volume_info(dh).get_corner_dart(i, 0)));

      for (int j=0; j<4; ++j)
        middledarts[i][j]=get_middle_edge(volume_info(dh).get_corner_dart(i, j));

      facedarts[i]=get_middle_face(middledarts[i][0]);
    }

    //display_corner_darts(m_lcc->template attribute<3>(dh)->get_corner_darts());
    //display_middle_darts(middledarts);
    //display_face_darts(facedarts);

    make_eight_hexahedron(volume_info(dh).get_corner_darts(),
                          middledarts, facedarts);
  }

  // Subdivide the hexahedron given by one of its dart.
  // If the given dart does not belong to an hexahedron, the lcc is not modified.
  void subdivide(Dart_descriptor dh)
  {
    // Logger logger("Linear_cell_complex_subdivision::subdivide");

    if(m_lcc->template attribute<3>(dh)==m_lcc->null_handle) return;

#ifdef LOG_DISTRIBUTED_LCC
    std::cout<<"[LCC_subdivision] : Split hexa (";
    std::vector<Point> eight_points;
    for (unsigned int i=0; i<8; ++i)
    {
      eight_points.push_back(m_lcc->point_of_vertex_attribute
                             (m_lcc->template attribute<3>(dh)->info().
                              get_corner_particle(*m_lcc, i)));
    }
    std::sort(eight_points.begin(), eight_points.end());
    for (unsigned int i=0; i<8; ++i)
    {
      std::cout<<eight_points[i];
      if (i<7) std::cout<<", ";
    }
    std::cout<<")"<<std::endl;
#endif

    // debug("Before to subdivide");

    subdivide_edges_of_hexahedron(dh);
    //assert(check_corner_darts());

    // debug("After subdivide edges");

    subdivide_faces_of_hexahedron(dh);
    //assert(check_corner_darts());

    // debug("After subdivide faces");

    subdivide_hexahedron(dh);
    assert(m_lcc->is_valid());
    //assert(check_corner_darts()); // TODO BUG ?

    // debug("After subdivide hexahedron");

#ifdef LOG_DISTRIBUTED_LCC
    std::cout<<"[LCC_subdivision] : Split hexa  END"<<std::endl;
#endif
  }

  /// Split a given edge: the edge is split in two.
  void split_edge_in_two(Dart_descriptor dh)
  {
    // Logger logger("Linear_cell_complex_subdivision::split_edge_in_two");

#ifdef LOG_DISTRIBUTED_LCC
    std::cout<<"[LCC_subdivision] : Split edge ("
            <<m_lcc->point(dh)<<", "
           <<m_lcc->point(m_lcc->opposite(dh))
          <<")";
    std::cout<<std::endl;
#endif
    typename LCC::size_type m=m_lcc->get_new_mark();
    m_lcc->negate_mark(m);

    m_lcc->template insert_barycenter_in_cell<1>(dh);

    m_lcc->negate_mark(m);
    // Now only new darts are marked

    for(typename LCC::template Dart_of_cell_range<1>::iterator
          it(m_lcc->template darts_of_cell<1>(dh).begin()),
          itend(m_lcc->template darts_of_cell<1>(dh).end());
        it!=itend; ++it)
    {
      if (m_lcc->is_marked(it, m))
      { // Here it is a new dart
        set_first_direction // set direction only of new darts
            (it, is_first_direction(m_lcc->template beta<0>(it)));

        int level=2*get_subdivision_level_of_edge(m_lcc->template beta<0>(it));
        set_subdivision_level_of_edge(it, level);
        set_subdivision_level_of_edge(m_lcc->template beta<0>(it), level);

        m_lcc->unmark(it, m);
      }
      else
      { // Here it is an old dart
        assert(m_lcc->is_marked(m_lcc->template beta<1>(it), m));

        set_first_direction // set direction only of new darts
            (m_lcc->template beta<1>(it), is_first_direction(it));

        int level=2*get_subdivision_level_of_edge(it);
        set_subdivision_level_of_edge(it, level);
        set_subdivision_level_of_edge(m_lcc->template beta<1>(it), level);

        m_lcc->unmark(m_lcc->template beta<1>(it), m);
      }
    }

    assert(m_lcc->is_whole_map_unmarked(m));
    m_lcc->free_mark(m);

    if (split_edge_in_two_functor()!=nullptr)
      split_edge_in_two_functor()(dh, m_lcc->template beta<1>(dh));
  }

protected:
  /// @return true iff dh is a corner dart.
  bool is_corner_dart(Dart_descriptor dh)
  { return m_lcc->is_marked(dh, m_corner_mark); }

  /// Mark dh as corner.
  void set_corner_dart(Dart_descriptor dh)
  { m_lcc->mark(dh, m_corner_mark); }

  /// Unmark dh as corner.
  void unset_corner_dart(Dart_descriptor dh)
  { m_lcc->unmark(dh, m_corner_mark); }

  /// @return true iff dh is a dart that belongs to an edge in the first direction.
  bool is_first_direction(Dart_descriptor dh)
  { return m_lcc->is_marked(dh, m_edge_first_direction); }

  /// Set dh as first direction.
  void set_first_direction(Dart_descriptor dh)
  { m_lcc->mark(dh, m_edge_first_direction); }

  /// Set dh as second direction.
  void set_second_direction(Dart_descriptor dh)
  { m_lcc->unmark(dh, m_edge_first_direction); }

  /// Mark dh as first direction or second depending on b.
  void set_first_direction(Dart_descriptor dh, bool b)
  { m_lcc->set_mark_to(dh, m_edge_first_direction, b); }

  /// @return true iff the edge containing dh is already subdivided.
  bool is_edge_already_subdivided(Dart_descriptor dh)
  {
    // An edge is already subdivided if beta1(dh) is not a corner dart.
    // @TODO if we want other subdivision scheme (split an edge in three...) we
    //     need to test if an edge is split in k by searching if there is a sum
    //     of edge size which is equal to 1/k of hierarchical edge size
    return !is_corner_dart(m_lcc->template beta<1>(dh));
  }

  /// @return true iff the face containing dh is already subdivided.
  bool is_face_already_subdivided(Dart_descriptor dh)
  {
    // A face is already subdivided if when following its border we change
    // of direction without passing through a corner dart.
    // @TODO same remark than for is_edge_already_subdivided if we want to
    // make different subdvision scheme, in this case we need to test if
    // the subdvision of the face is "compatible".

    Dart_descriptor cur = dh;
    bool prevdir=is_first_direction(dh);
    do
    {
      cur=m_lcc->template beta<1>(cur);
      if (is_first_direction(cur)!=prevdir &&
          !is_corner_dart(cur)) return true;
      prevdir=is_first_direction(cur);
    }
    while (cur!=dh);
    return false;
  }

  /** @return the subdivision level of the edge containing dart dh.
   *  Each dart is associated with an level_size (unsigned int) giving the
   *     subdivision level l of its edge.
   *  This level is always regarding the hierarchical edge, ie the boundary of the hexahedra
   *  that contains this edge.
   *  When an edge is cut in k, we multiply its subdivision level by k.
   *  When a hexadra is cut in k, we divide the subdivision level of its edge by k
   *   (to take into account the fact that the new hierarchical edge is now shorter).
   */
  level_size get_subdivision_level_of_edge(Dart_descriptor dh)
  {
    if (!m_dart_levels.is_defined(dh) ) return 1;
    return m_dart_levels[dh];
  }

  /// Set the subdivision level of edge for dart dh to l.
  ///   (@see get_subdivision_level_of_edge)
  void set_subdivision_level_of_edge(Dart_descriptor dh, level_size l)
  { m_dart_levels[dh] = l; }

  /// Multiply the subdivision level of edge for dart dh by l.
  ///   (@see get_subdivision_level_of_edge)
  void multiply_subdivision_level_of_edge(Dart_descriptor dh, level_size l)
  {
    level_size& level=m_dart_levels[dh];
    level *= l;
  }

  /// Divide the subdivision level of edge for dart dh by l
  ///   (@see get_subdivision_level_of_edge)
  void divide_subdivision_level_of_edge(Dart_descriptor dh, level_size l)
  {
    level_size& level=m_dart_levels[dh];
    assert(level % l==0);
    level /= l;
  }

  /// @return the next edge dart after dh by following the same direction.
  /// here we suppose edges traversing the followed edge are labeled with another direction.
  Dart_descriptor get_next_edge_dart(Dart_descriptor dh)
  {
    Dart_descriptor res = m_lcc->template beta<1>(dh);

    assert(!is_corner_dart(res));

    // We are sure that the following loop will finish by construction
    // of border darts.
    // Moreover with the current version (hexa only one subdvision), we are
    // sure we can make only one loop.
    while (is_first_direction(res)!=is_first_direction(dh))
    {
      res = m_lcc->template beta<2,1>(res);
      assert(res!=m_lcc->null_dart_handle);
      assert(res!=m_lcc->template beta<1>(dh));
    }
    return res;
  }

  /// @return the previous border dart before dh.
  Dart_descriptor get_previous_edge_dart(Dart_descriptor dh)
  {
    assert(!is_corner_dart(dh));

    Dart_descriptor res = m_lcc->template beta<0>(dh);

    // We are sure that the following loop will finish by construction
    // of border darts.
    while (is_first_direction(res)!=is_first_direction(dh))
    {
      res = m_lcc->template beta<2,0>(res);
      assert(res!=m_lcc->null_dart_handle);
      assert(res!=m_lcc->template beta<0>(dh));
    }
    return res;
  }

  /// @return the next corner dart after dh.
  Dart_descriptor get_next_corner_dart(Dart_descriptor dh)
  {
    assert(is_corner_dart(dh));
    while (!is_corner_dart(m_lcc->template beta<1>(dh)))
    {
      dh = get_next_edge_dart(dh);
    }
    return m_lcc->template beta<1>(dh);
  }

  // Compute the subdvision level of hierarchical edge but for opposite
  // hexahedra (ie the one containing beta3(dh).
  // This function is used before to split a face in 4. In this case, we know
  // that there is only one opposite hexahedra (otherwise the face would be
  // already subdivided). Moreover the current face has 4 corner darts.
  level_size compute_subdivision_level_opposite_edge(Dart_descriptor dh)
  {
    assert(!m_lcc->template is_free<3>(dh));
    assert(is_corner_dart(dh));

    CFraction sum(0, 1);
    do
    {
      sum += get_subdivision_level_of_edge(m_lcc->template beta<3>(dh));
      dh = m_lcc->template beta<1>(dh);
    }
    while (!is_corner_dart(dh));

    assert(sum.get_num()==1);
    return sum.get_den();
  }

  /// @return the edge in the middle of the hierarchical edge containing dart dh.
  /// @pre dh is a corner dart;
  ///      the edge containing dh is already subdivided.
  Dart_descriptor get_middle_edge(Dart_descriptor dh)
  {
    assert(is_corner_dart(dh));

    // Half of the subdivision level of the hierarchical edge
    CFraction size_he(1, 2);

    // Fraction for the sum of edge sizes
    CFraction sum(0, 1);
    while(sum<size_he)
    {
      sum += get_subdivision_level_of_edge(dh);
      dh = get_next_edge_dart(dh);
    }

    assert(sum==size_he);
    return dh;
  }

  /** @return the edge in the middle of a hierarchical edge containing dart dh
   *          given the subdvision level of the hierarchical edge.
   *  @pre the edge containing dh is already subdivided.
   */
  Dart_descriptor get_middle_edge(Dart_descriptor dh, level_size level)
  {
    // Half of the subdivision level of the hierarchical edge
    CFraction size_he(1, 2*level);

    // Fraction for the sum of edge sizes
    CFraction sum(0, 1);
    while(sum<size_he)
    {
      sum += get_subdivision_level_of_edge(dh);
      dh = m_lcc->template beta<1>(dh);
    }

    assert(sum==size_he);
    return dh;
  }

  /// @return the edge in the middle of the hierarchical face containing dart dh.
  /// @pre dh is a border dart in the middle of one hierarchical edge;
  ///      the face containing dh is already subdivided.
  Dart_descriptor get_middle_face(Dart_descriptor dh)
  {
    // Half of the subdivision level of the hierarchical edge
    CFraction size_he(1, 2);

    // Fraction for the sum of edge sizes
    CFraction sum(0, 1);
    dh = m_lcc->template beta<0,2>(dh);

    while(sum<size_he)
    {
      sum += get_subdivision_level_of_edge(dh);
      dh = get_next_edge_dart(dh);
    }

    assert(sum==size_he);
    return dh;
  }

  /// Split a given face in four.
  /// @pre: the face must not be already subdivided;
  /// @pre: the edges of the face must be already subdivided.
  /// @pre: dh1,...,dh4 are the four corner darts of the face
  void split_face_in_four(Dart_descriptor dh1,
                          Dart_descriptor dh2,
                          Dart_descriptor dh3,
                          Dart_descriptor dh4)
  {
#ifdef LOG_DISTRIBUTED_LCC
    std::cout<<"[LCC_subdivision] : Split face ("
            <<m_lcc->point(dh1)<<", "
           <<m_lcc->point(dh2)<<", "
          <<m_lcc->point(dh3)<<", "
         <<m_lcc->point(dh4)<<")"<<std::endl;
#endif

    assert(!is_face_already_subdivided(dh1));

    assert(is_corner_dart(dh1));
    assert(is_corner_dart(dh2));
    assert(is_corner_dart(dh3));
    assert(is_corner_dart(dh4));

    assert(is_first_direction(dh1)==is_first_direction(dh3));
    assert(is_first_direction(dh2)==is_first_direction(dh4));

    level_size oppositelevel=0;
    if (!m_lcc->template is_free<3>(dh1))
    {
      oppositelevel=2*compute_subdivision_level_opposite_edge(dh1);
      assert(oppositelevel==2*compute_subdivision_level_opposite_edge(dh2));
      assert(oppositelevel==2*compute_subdivision_level_opposite_edge(dh3));
      assert(oppositelevel==2*compute_subdivision_level_opposite_edge(dh4));
    }

    Dart_descriptor mh1 = get_middle_edge(dh1);
    Dart_descriptor mh2 = get_middle_edge(dh2);
    Dart_descriptor mh3 = get_middle_edge(dh3);
    Dart_descriptor mh4 = get_middle_edge(dh4);

    // Here all the new created dart are not border dart nor corner dart.

    Dart_descriptor e1 = m_lcc->insert_cell_1_in_cell_2(mh1, mh3);
    m_lcc->template insert_barycenter_in_cell<1>(e1);

    set_subdivision_level_of_edge(e1,                            2);
    set_subdivision_level_of_edge(m_lcc->template beta<2>  (e1), 2);
    set_subdivision_level_of_edge(m_lcc->template beta<1>  (e1), 2);
    set_subdivision_level_of_edge(m_lcc->template beta<1,2>(e1), 2);

    set_first_direction(e1,                            is_first_direction(dh2));
    set_first_direction(m_lcc->template beta<1>  (e1), is_first_direction(dh2));
    set_first_direction(m_lcc->template beta<2>  (e1), is_first_direction(dh2));
    set_first_direction(m_lcc->template beta<1,2>(e1), is_first_direction(dh2));

    if (!m_lcc->template is_free<3>(e1))
    {
      set_subdivision_level_of_edge(m_lcc->template beta<3>    (e1), oppositelevel);
      set_subdivision_level_of_edge(m_lcc->template beta<2,3>  (e1), oppositelevel);
      set_subdivision_level_of_edge(m_lcc->template beta<1,3>  (e1), oppositelevel);
      set_subdivision_level_of_edge(m_lcc->template beta<1,2,3>(e1), oppositelevel);

      set_first_direction(m_lcc->template beta<3>    (e1),
                          is_first_direction(m_lcc->template beta<3>(dh2)));
      set_first_direction(m_lcc->template beta<1,3>  (e1),
                          is_first_direction(m_lcc->template beta<3>(dh2)));
      set_first_direction(m_lcc->template beta<2,3>  (e1),
                          is_first_direction(m_lcc->template beta<3>(dh2)));
      set_first_direction(m_lcc->template beta<1,2,3>(e1),
                          is_first_direction(m_lcc->template beta<3>(dh2)));
    }

    Dart_descriptor e2 =
        m_lcc->insert_cell_1_in_cell_2(m_lcc->template beta<1>(e1), mh2);
    Dart_descriptor e3 =
        m_lcc->insert_cell_1_in_cell_2(m_lcc->template beta<2>(e1), mh4);

    set_subdivision_level_of_edge(e2,                          2);
    set_subdivision_level_of_edge(m_lcc->template beta<2>(e2), 2);
    set_subdivision_level_of_edge(e3,                          2);
    set_subdivision_level_of_edge(m_lcc->template beta<2>(e3), 2);

    set_first_direction(e2                         , is_first_direction(dh1));
    set_first_direction(m_lcc->template beta<2>(e2), is_first_direction(dh1));
    set_first_direction(e3                         , is_first_direction(dh1));
    set_first_direction(m_lcc->template beta<2>(e3), is_first_direction(dh1));

    if (!m_lcc->template is_free<3>(e1))
    {
      set_subdivision_level_of_edge(m_lcc->template beta<3>(e2),   oppositelevel);
      set_subdivision_level_of_edge(m_lcc->template beta<2,3>(e2), oppositelevel);
      set_subdivision_level_of_edge(m_lcc->template beta<3>(e3),   oppositelevel);
      set_subdivision_level_of_edge(m_lcc->template beta<2,3>(e3), oppositelevel);

      set_first_direction(m_lcc->template beta<3>(e2),
                          is_first_direction(m_lcc->template beta<3>(dh1)));
      set_first_direction(m_lcc->template beta<2,3>(e2),
                          is_first_direction(m_lcc->template beta<3>(dh1)));
      set_first_direction(m_lcc->template beta<3>(e3),
                          is_first_direction(m_lcc->template beta<3>(dh1)));
      set_first_direction(m_lcc->template beta<2,3>(e3),
                          is_first_direction(m_lcc->template beta<3>(dh1)));
    }

    if (split_face_in_four_functor()!=nullptr)
      split_face_in_four_functor()(e1, e2, e3, m_lcc->template beta<1,2>(e2));
  }

  /** Split a given face in four, for a face which does not correspond to the
   *  border of an hexahedron.
   *  Note that the edges of the face are not necessarily subdivided.
   */
  void split_face_in_four(Dart_descriptor dh1)
  {
    // Logger logger("Linear_cell_complex_subdivision::split_face_in_four");

    level_size size1a, size2a, size3a, size4a;
    level_size size1b, size2b, size3b, size4b;

    if (is_first_direction(dh1)!=is_first_direction(m_lcc->template beta<1>(dh1)))
      split_edge_in_two(dh1);
    Dart_descriptor dh2 = compute_length_of_face_border(dh1, size1a, size1b);

    if (is_first_direction(dh2)!=is_first_direction(m_lcc->template beta<1>(dh2)))
      split_edge_in_two(dh2);
    Dart_descriptor dh3 = compute_length_of_face_border(dh2, size2a, size2b);

    if (is_first_direction(dh3)!=is_first_direction(m_lcc->template beta<1>(dh3)))
      split_edge_in_two(dh3);
    Dart_descriptor dh4 = compute_length_of_face_border(dh3, size3a, size3b);

    if (is_first_direction(dh4)!=is_first_direction(m_lcc->template beta<1>(dh4)))
      split_edge_in_two(dh4);
    compute_length_of_face_border(dh4, size4a, size4b);

    assert(size1a==size3a);
    assert(size2a==size4a);
    assert(size1b==size3b);
    assert(size2b==size4b);

    assert(is_first_direction(dh1)==is_first_direction(dh3));
    assert(is_first_direction(dh2)==is_first_direction(dh4));

    Dart_descriptor mh1 = get_middle_edge(dh1, size1a);
    Dart_descriptor mh2 = get_middle_edge(dh2, size2a);
    Dart_descriptor mh3 = get_middle_edge(dh3, size3a);
    Dart_descriptor mh4 = get_middle_edge(dh4, size4a);

#ifdef LOG_DISTRIBUTED_LCC
    std::cout<<"[LCC_subdivision] : Split face(b) ("
            <<m_lcc->point(dh1)<<", "
           <<m_lcc->point(dh2)<<", "
          <<m_lcc->point(dh3)<<", "
         <<m_lcc->point(dh4)<<")";
    std::cout<<std::endl;
#endif

    // Here all the new created dart are not border dart nor corner dart.

    Dart_descriptor e1 = m_lcc->insert_cell_1_in_cell_2(mh1, mh3);
    m_lcc->template insert_barycenter_in_cell<1>(e1);

    set_subdivision_level_of_edge(e1,                            2*size2a);
    set_subdivision_level_of_edge(m_lcc->template beta<2>  (e1), 2*size2a);
    set_subdivision_level_of_edge(m_lcc->template beta<1>  (e1), 2*size2a);
    set_subdivision_level_of_edge(m_lcc->template beta<1,2>(e1), 2*size2a);

    set_first_direction(e1,                            is_first_direction(dh2));
    set_first_direction(m_lcc->template beta<1>  (e1), is_first_direction(dh2));
    set_first_direction(m_lcc->template beta<2>  (e1), is_first_direction(dh2));
    set_first_direction(m_lcc->template beta<1,2>(e1), is_first_direction(dh2));

    if (!m_lcc->template is_free<3>(e1))
    {
      set_subdivision_level_of_edge(m_lcc->template beta<3>    (e1), 2*size2b);
      set_subdivision_level_of_edge(m_lcc->template beta<2,3>  (e1), 2*size2b);
      set_subdivision_level_of_edge(m_lcc->template beta<1,3>  (e1), 2*size2b);
      set_subdivision_level_of_edge(m_lcc->template beta<1,2,3>(e1), 2*size2b);

      set_first_direction(m_lcc->template beta<3>    (e1),
                          is_first_direction(m_lcc->template beta<3>(dh2)));
      set_first_direction(m_lcc->template beta<1,3>  (e1),
                          is_first_direction(m_lcc->template beta<3>(dh2)));
      set_first_direction(m_lcc->template beta<2,3>  (e1),
                          is_first_direction(m_lcc->template beta<3>(dh2)));
      set_first_direction(m_lcc->template beta<1,2,3>(e1),
                          is_first_direction(m_lcc->template beta<3>(dh2)));
    }

    Dart_descriptor e2 =
        m_lcc->insert_cell_1_in_cell_2(m_lcc->template beta<1>(e1), mh2);
    Dart_descriptor e3 =
        m_lcc->insert_cell_1_in_cell_2(m_lcc->template beta<2>(e1), mh4);

    set_subdivision_level_of_edge(e2,                          2*size1a);
    set_subdivision_level_of_edge(m_lcc->template beta<2>(e2), 2*size1a);
    set_subdivision_level_of_edge(e3,                          2*size1a);
    set_subdivision_level_of_edge(m_lcc->template beta<2>(e3), 2*size1a);

    set_first_direction(e2                         , is_first_direction(dh1));
    set_first_direction(m_lcc->template beta<2>(e2), is_first_direction(dh1));
    set_first_direction(e3                         , is_first_direction(dh1));
    set_first_direction(m_lcc->template beta<2>(e3), is_first_direction(dh1));

    if (!m_lcc->template is_free<3>(e1))
    {
      set_subdivision_level_of_edge(m_lcc->template beta<3>(e2),   2*size1b);
      set_subdivision_level_of_edge(m_lcc->template beta<2,3>(e2), 2*size1b);
      set_subdivision_level_of_edge(m_lcc->template beta<3>(e3),   2*size1b);
      set_subdivision_level_of_edge(m_lcc->template beta<2,3>(e3), 2*size1b);

      set_first_direction(m_lcc->template beta<3>(e2),
                          is_first_direction(m_lcc->template beta<3>(dh1)));
      set_first_direction(m_lcc->template beta<2,3>(e2),
                          is_first_direction(m_lcc->template beta<3>(dh1)));
      set_first_direction(m_lcc->template beta<3>(e3),
                          is_first_direction(m_lcc->template beta<3>(dh1)));
      set_first_direction(m_lcc->template beta<2,3>(e3),
                          is_first_direction(m_lcc->template beta<3>(dh1)));
    }

    if (split_face_in_four_functor()!=nullptr)
      split_face_in_four_functor()(e1, e2, e3, m_lcc->template beta<1,2>(e2));
  }

  /** Identify the two border d1 and d2:
   *  d1 is a dart from the border of the hexahedra (thus the edge could be subdivided)
   *  d2 is a dart from the new created pattern (thus the edge is never subdivided).
   *  2-link (beta2(d1),beta3(d2)) and (d1,d2); while subdividing edge d2 as necessary
   */
  void identify(Dart_descriptor d1, Dart_descriptor d2)
  {
    set_corner_dart(d1);
    set_corner_dart(m_lcc->template beta<2,1>(d1));

    bool cont = true;
    do
    {
      // std::cout<<m_lcc->point(m_lcc->other_extremity(d1))<<" et "<<m_lcc->point(d2)<<std::endl;
      if (m_lcc->vertex_attribute(m_lcc->other_extremity(d1))!=m_lcc->vertex_attribute(d2))
      {
        // We insert a vertex in d2, sharing its vertex attribute with the second extremity of d1
        m_lcc->insert_cell_0_in_cell_1(d2, m_lcc->vertex_attribute(m_lcc->other_extremity(d1)));
        d2=m_lcc->template beta<1>(d2);

        // In this configuration, we know that edge containing d2 has only two darts: d2 and b3(d2)
        // new darts (after insertion and affectation of d2) are thus d2 and b03(d2)
        set_first_direction(d2,
                            is_first_direction(m_lcc->template beta<0>(d2)));
        set_first_direction(m_lcc->template beta<0,3>(d2),
                            is_first_direction(m_lcc->template beta<3>(d2)));

        cont=true;
      }
      else
      {
        set_corner_dart(m_lcc->template beta<2>(d1));

        cont = false;
      }

      assert(get_subdivision_level_of_edge(m_lcc->template beta<2>(d1))==
             get_subdivision_level_of_edge(d1));

      m_lcc->template link_beta_for_involution<2>(m_lcc->template beta<2>(d1),
                                                  m_lcc->template beta<3>(d2));
      m_lcc->template link_beta_for_involution<2>(d1, d2);

    // Replace the 4 following lines by the TODO commented ones (if correct ???)
      m_lcc->template set_dart_attribute<3>(m_lcc->template beta<3>(d2), m_lcc->null_handle);
      m_lcc->template set_dart_attribute<3>(d2, m_lcc->null_handle);

      set_subdivision_level_of_edge(m_lcc->template beta<3>(d2),
                                    get_subdivision_level_of_edge(d1)/2);

      set_subdivision_level_of_edge(d2,
                                    get_subdivision_level_of_edge(d1)/2);

      /* TODO This should work if we comment line 980
       *  // if (m_lcc->template attribute<3>(itvol)!=m_lcc->null_handle)
       * but it does not.
      set_subdivision_level_of_edge(m_lcc->template beta<3>(d2),
                                    get_subdivision_level_of_edge(d1));

      set_subdivision_level_of_edge(d2,
                                    get_subdivision_level_of_edge(d1));
      */

      if (cont)
      {
        d1=get_next_edge_dart(d1);
        d2=m_lcc->template beta<0>(d2);
      }
    }
    while (cont);
  }

  void mark_hexahedron(Dart_descriptor dh)
  {
    set_first_direction(dh);
    set_first_direction(m_lcc->template beta<1, 1>(dh));

    set_first_direction(m_lcc->template beta<2>(dh));
    set_first_direction(m_lcc->template beta<2, 1, 1>(dh));

    set_first_direction(m_lcc->template beta<1, 1, 2>(dh));
    set_first_direction(m_lcc->template beta<1, 1, 2, 1, 1>(dh));

    set_first_direction(m_lcc->template beta<2, 1, 1, 2>(dh));
    set_first_direction(m_lcc->template beta<2, 1, 1, 2, 1, 1>(dh));

    set_first_direction(m_lcc->template beta<1, 2, 1>(dh));
    set_first_direction(m_lcc->template beta<1, 2, 0>(dh));

    set_first_direction(m_lcc->template beta<0, 2, 1>(dh));
    set_first_direction(m_lcc->template beta<0, 2, 0>(dh));
  }

  void initialize_face_corners(Dart_descriptor dh,
                               typename LCC::template Attribute_handle<3>::type a,
                               unsigned int num_face)
  {
    assert(num_face<6);
    for (unsigned int i=0; i<4; i++)
    {
      volume_info(a).set_corner_dart(num_face, i, dh);
      dh=get_next_corner_dart(dh);
    }
  }

  /** Update the 3-attribute of all dart of 3-cell(dh) to a; and divide by two
   *  all the subdivision levels of the old darts (indeed the new ones are
   *  initialized to 1; we know that a dart is old because its 3-attribute
   *  is m_lcc->null_handle). Then recompute all corner darts of the new hexahedron.
   */
  void update_new_hexahedron(Dart_descriptor dh,
                             typename LCC::template Attribute_handle<3>::type a,
                             unsigned char newlevel)
  {
    m_volume_info.resize(m_lcc->template attributes<3>().capacity());
    volume_info(a).set_subdivision_level(newlevel);

    for (typename LCC::template Dart_of_cell_range<3>::iterator
           itvol   =m_lcc->template darts_of_cell<3>(dh).begin(),
           itvolend=m_lcc->template darts_of_cell<3>(dh).end();
         itvol!=itvolend ; ++itvol)
    {
      if (m_lcc->template attribute<3>(itvol)!=m_lcc->null_handle)
      { divide_subdivision_level_of_edge(itvol, 2); }

      m_lcc->template set_dart_attribute<3>(itvol, a);
    }

    // Dart corners
    for (unsigned int i=0; i<4; i++)
    {
      initialize_face_corners(dh, a, i);

      // Move to an other face
      dh=m_lcc->template beta<0,2>(volume_info(a).get_corner_dart(i, 3));
    }//4 faces

    // Move to the down face
    dh=m_lcc->template beta<0,2>(volume_info(a).get_corner_dart(0, 0));
    initialize_face_corners(dh, a, 4);

    // Move to the up face
    dh=m_lcc->template beta<0,2>(volume_info(a).get_corner_dart(0, 2));
    initialize_face_corners(dh, a, 5);
  }

  Dart_descriptor make_quadrangle(Vertex_attribute_handle h0,
                              Vertex_attribute_handle h1,
                              Vertex_attribute_handle h2,
                              Vertex_attribute_handle h3)
  {
    Dart_descriptor d1 = m_lcc->make_combinatorial_polygon(4);

    m_lcc->set_vertex_attribute_of_dart(d1,h0);
    m_lcc->set_vertex_attribute_of_dart(m_lcc->template beta<1>(d1), h1);
    m_lcc->set_vertex_attribute_of_dart(m_lcc->template beta<1, 1>(d1), h2);
    m_lcc->set_vertex_attribute_of_dart(m_lcc->template beta<0>(d1), h3);

    set_first_direction(d1);
    set_first_direction(m_lcc->template beta<1,1>(d1));

    return d1;
  }

  Dart_descriptor sew_combinatorial_hexahedron(Dart_descriptor d1,
                                           Dart_descriptor d2,
                                           Dart_descriptor d3,
                                           Dart_descriptor d4,
                                           Dart_descriptor d5,
                                           Dart_descriptor d6)
  {
    if (d1!=m_lcc->null_handle)
    {
      if (d4!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(d1,
                                              m_lcc->beta(d4, 1, 1), 2);
      if (d6!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(m_lcc->beta(d1, 1),
                                              m_lcc->beta(d6, 0)    , 2);
      if (d2!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(m_lcc->beta(d1, 1, 1),
                                              d2                  , 2);
      if (d5!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(m_lcc->beta(d1, 0),
                                              d5                  , 2);
    }

    if (d3!=m_lcc->null_handle)
    {
      if (d2!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(d3,
                                              m_lcc->beta(d2, 1, 1), 2);
      if (d6!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(m_lcc->beta(d3, 1),
                                              m_lcc->beta(d6, 1)         , 2);
      if (d4!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(m_lcc->beta(d3, 1, 1),
                                              d4                  , 2);
      if (d5!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(m_lcc->beta(d3, 0),
                                              m_lcc->beta(d5, 1, 1), 2);
    }

    if (d6!=m_lcc->null_handle)
    {
      if (d4!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(d6,
                                              m_lcc->beta(d4, 1)         , 2);
      if (d2!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(m_lcc->beta(d6, 1, 1),
                                              m_lcc->beta(d2, 1)         , 2);
    }

    if (d5!=m_lcc->null_handle)
    {
      if (d4!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(m_lcc->beta(d5, 0),
                                              m_lcc->beta(d4, 0)         , 2);
      if (d2!=m_lcc->null_handle)
        m_lcc->basic_link_beta_for_involution(m_lcc->beta(d5, 1),
                                              m_lcc->beta(d2, 0)         , 2);
    }

    if (d1!=m_lcc->null_handle)
      return d1;

    if (d4!=m_lcc->null_handle)
      return m_lcc->template beta<1,1>(d4);

    assert(d2!=m_lcc->null_handle);
    return d2;
  }

  Dart_descriptor make_partial_hexahedron(Vertex_attribute_handle h0,
                                      Vertex_attribute_handle h1,
                                      Vertex_attribute_handle h2,
                                      Vertex_attribute_handle h3,
                                      Vertex_attribute_handle h4,
                                      Vertex_attribute_handle h5,
                                      Vertex_attribute_handle h6,
                                      Vertex_attribute_handle h7,
                                      bool withf0,
                                      bool withf1,
                                      bool withf2,
                                      bool withf3,
                                      bool withf4,
                                      bool withf5)
  {
    Dart_descriptor d0=m_lcc->null_handle, d1=m_lcc->null_handle;
    Dart_descriptor d2=m_lcc->null_handle, d3=m_lcc->null_handle;
    Dart_descriptor d4=m_lcc->null_handle, d5=m_lcc->null_handle;
    if(withf0) { d0=make_quadrangle(h0, h5, h6, h1); }
    if(withf1) { d1=make_quadrangle(h1, h6, h7, h2); }
    if(withf2) { d2=make_quadrangle(h2, h7, h4, h3); }
    if(withf3) { d3=make_quadrangle(h3, h4, h5, h0); }
    if(withf4) { d4=make_quadrangle(h0, h1, h2, h3); }
    if(withf5) { d5=make_quadrangle(h5, h4, h7, h6); }

    return sew_combinatorial_hexahedron(d0, d1, d2, d3, d4, d5);
  }

  /** Create a 2x2x2 hexahedral grid given 8 points.
   * \verbatim
   *
   * Vertices of all hexahedron:
   *         y----z----aa
   *        /|   /|   /|
   *       v----w----x |
   *      /| p-/--q-/|-r
   *     s----t----u |/|
   *     | m--|-n--|-o |
   *     |/| g|---h|/|-i
   *     j----k----l |/
   *     | d--|-e--|-f
   *     |/   |/   |/
   *     a----b----c
   *
   * Number of each hexahedron:
   *       6----7
   *      /|   /|
   *     4----5 |
   *     | 2--|-3
   *     |/   |/
   *     0----1
   *
   */
  void make_eight_hexahedron(Dart_descriptor (&cd) [6][4], // cd=array of corner darts
                             Dart_descriptor (&md) [6][4], // md=array of middle darts
                             Dart_descriptor (&fd) [6]) // fd=array of middle darts
  {
    // Logger logger("Linear_cell_complex_subdivision::make_eight_hexahedron");

    // We negate m_corner_mark so that all new darts will be marked as corners.
    m_lcc->negate_mark(m_corner_mark);

    Vertex_attribute_handle center=m_lcc->create_vertex_attribute
      ( typename LCC::Traits::Construct_midpoint()(m_lcc->point(fd[0]),
                                                   m_lcc->point(fd[2])));

#define VH(d) m_lcc->vertex_attribute(d)

    // 1) Create 8 partial cubes
    Dart_descriptor d0 =
        make_partial_hexahedron(VH(cd[4][0]), VH(md[4][0]),
                                VH(fd[4]), VH(md[4][3]),
                                VH(fd[3]), VH(md[0][0]),
                                VH(fd[0]), center,
                                false, true, true, false, false, true );
    Dart_descriptor d1 =
        make_partial_hexahedron(VH(md[4][0]), VH(cd[4][1]),
                                VH(md[4][1]), VH(fd[4]),
                                center, VH(fd[0]),
                                VH(md[0][2]), VH(fd[1]),
                                false, false, true, true, false, true);
    Dart_descriptor d2 =
        make_partial_hexahedron(VH(md[4][3]), VH(fd[4]),
                                VH(md[4][2]), VH(cd[4][3]),
                                VH(md[2][2]), VH(fd[3]),
                                center, VH(fd[2]),
                                true, true, false, false, false, true);
    Dart_descriptor d3 =
        make_partial_hexahedron(VH(fd[4]), VH(md[4][1]),
                                VH(cd[4][2]), VH(md[4][2]),
                                VH(fd[2]), center,
                                VH(fd[1]), VH(md[1][2]),
                                true, false, false, true, false, true);
    Dart_descriptor d4 =
        make_partial_hexahedron(VH(md[0][0]), VH(fd[0]),
                                center, VH(fd[3]),
                                VH(md[5][1]), VH(cd[5][1]),
                                VH(md[5][0]), VH(fd[5]),
                                false, true, true, false, true, false);
    Dart_descriptor d5 =
        make_partial_hexahedron(VH(fd[0]), VH(md[0][2]),
                                VH(fd[1]), center,
                                VH(fd[5]), VH(md[5][0]),
                                VH(cd[5][0]), VH(md[5][3]),
                                false, false, true, true, true, false);
    Dart_descriptor d6 =
        make_partial_hexahedron(VH(fd[3]), center,
                                VH(fd[2]), VH(md[2][2]),
                                VH(cd[5][2]), VH(md[5][1]),
                                VH(fd[5]), VH(md[5][2]),
                                true, true, false, false, true, false);
    Dart_descriptor d7 =
        make_partial_hexahedron(center, VH(fd[1]),
                                VH(md[1][2]), VH(fd[2]),
                                VH(md[5][2]), VH(fd[5]),
                                VH(md[5][3]), VH(cd[5][3]),
                                true, false, false, true, true, false);

#undef VH

    // We re-negate corner mark.
    m_lcc->negate_mark(m_corner_mark);

    // 2) 3-sew the cubes between them
    // right face of h0 with left face of h1
    m_lcc->template topo_sew_for_involution<3>(d0, d1);
    // right face of h2 with left face of h3
    m_lcc->template topo_sew_for_involution<3>(m_lcc->template beta<1,1,2>(d2), m_lcc->template beta<2>(d3));
    // right face of h4 with left face of h5
    m_lcc->template topo_sew_for_involution<3>(d4, d5);
    // right face of h6 with left face of h7
    m_lcc->template topo_sew_for_involution<3>(m_lcc->template beta<1,1,2>(d6),m_lcc->template beta<2>(d7));

    // behind face of h0 with front face of h2
    m_lcc->template topo_sew_for_involution<3>(m_lcc->template beta<1,1,2,1,1>(d0), d2);
    // behind face of h1 with front face of h3
    m_lcc->template topo_sew_for_involution<3>(m_lcc->template beta<1,1,2>(d1), d3);
    // behind face of h4 with front face of h6
    m_lcc->template topo_sew_for_involution<3>(m_lcc->template beta<1,1,2,1,1>(d4), d6);
    // behind face of h5 with front face of h7
    m_lcc->template topo_sew_for_involution<3>(m_lcc->template beta<1,1,2>(d5), d7);

    // top face of h0 with bottom face of h4
    m_lcc->template topo_sew_for_involution<3>(m_lcc->template beta<1,2>(d0), m_lcc->template beta<0,2>(d4));
    // top face of h1 with bottom face of h5
    m_lcc->template topo_sew_for_involution<3>(m_lcc->template beta<0,2>(d1), m_lcc->template beta<1,2>(d5));
    // top face of h2 with bottom face of h6
    m_lcc->template topo_sew_for_involution<3>(m_lcc->template beta<1,2>(d2), m_lcc->template beta<0,2>(d6));
    // top face of h3 with bottom face of h7
    m_lcc->template topo_sew_for_involution<3>(m_lcc->template beta<1,2>(d3), m_lcc->template beta<0,2>(d7));

    // Keep one dart per inner new face
    Dart_descriptor nfd[12];

    // 4 right faces of the left volumes
    nfd[0]=m_lcc->template beta<1,1>(d0);
    nfd[1]=m_lcc->template beta<0>(d4);
    nfd[2]=m_lcc->template beta<1,1,2>(d6);
    nfd[3]=m_lcc->template beta<1,1,2,1>(d2);

    // 4 top faces of the low level volumes
    nfd[4]=m_lcc->template beta<1,2>(d0);
    nfd[5]=m_lcc->template beta<1,2>(d2);
    nfd[6]=m_lcc->template beta<1,2,1>(d3);
    nfd[7]=m_lcc->template beta<0,2,1>(d1);

    // 4 behind faces of the front volumes
    nfd[8]=m_lcc->template beta<1,1,2,1>(d0);
    nfd[9]=m_lcc->template beta<1,1,2>(d4);
    nfd[10]=m_lcc->template beta<1,1,2,1>(d5);
    nfd[11]=m_lcc->template beta<1,1,2>(d1);

    // Identification of cells: all created darts will be marked border darts
    Dart_descriptor toidentify[6][4];

    toidentify[0][0]=m_lcc->template beta<1,2,1>(d0);
    toidentify[0][1]=d4;
    toidentify[0][2]=m_lcc->template beta<1,2,1>(d5);
    toidentify[0][3]=d1;

    toidentify[1][0]=m_lcc->template beta<0,2,1,1>(d1);
    toidentify[1][1]=m_lcc->template beta<1,1,2,1,1>(d5);
    toidentify[1][2]=m_lcc->template beta<0,2,1>(d7);
    toidentify[1][3]=m_lcc->template beta<1,1>(d3);

    toidentify[2][0]=m_lcc->template beta<1,2,1,1>(d3);
    toidentify[2][1]=m_lcc->template beta<2,1,1>(d7);
    toidentify[2][2]=m_lcc->template beta<0,2,1,1>(d6);
    toidentify[2][3]=m_lcc->template beta<1,1,2,1,1>(d2);

    toidentify[3][0]=m_lcc->template beta<1,2,1>(d2);
    toidentify[3][1]=d6;
    toidentify[3][2]=m_lcc->template beta<0,2,1,1>(d4);
    toidentify[3][3]=m_lcc->template beta<1,1,2,1,1>(d0);

    toidentify[4][0]=m_lcc->template beta<0>(d0);
    toidentify[4][1]=m_lcc->template beta<1,1,2,1>(d1);
    toidentify[4][2]=m_lcc->template beta<2,1>(d3);
    toidentify[4][3]=m_lcc->template beta<0>(d2);

    toidentify[5][0]=m_lcc->template beta<0>(d5);
    toidentify[5][1]=m_lcc->template beta<1,1,2,1>(d4);
    toidentify[5][2]=m_lcc->template beta<1,1,2,1>(d6);
    toidentify[5][3]=m_lcc->template beta<1>(d7);

    for (int i=0; i<6; ++i)
      for (int j=0; j<4; ++j)
        identify(m_lcc->template beta<0,2>(md[i][j]), toidentify[i][j]);

    unsigned char nextlevel=volume_info(cd[0][0]).get_subdvision_level()+1;

    // Update new hexahedra: perhaps the updating of corner darts can be
    // optimized because we know some corners during the subdvision steps
    // (but more complex code).
    Dart_descriptor cd00=cd[0][0]; // Store this dart because the reserve can change
    // the array of corner darts, and thus invalidate cd.
    m_volume_info.reserve(m_volume_info.size()+7);
    update_new_hexahedron(cd00, m_lcc->template attribute<3>(cd00), nextlevel); // for hexa h0
    update_new_hexahedron(m_lcc->template beta<3,2>(d0), m_lcc->template create_attribute<3>(), nextlevel); // for hexa h1
    update_new_hexahedron(d2, m_lcc->template create_attribute<3>(), nextlevel); // for hexa h2
    update_new_hexahedron(d3, m_lcc->template create_attribute<3>(), nextlevel); // for hexa h3
    update_new_hexahedron(md[0][0], m_lcc->template create_attribute<3>(), nextlevel); // for hexa h4
    update_new_hexahedron(m_lcc->template beta<3,2>(d4), m_lcc->template create_attribute<3>(), nextlevel); // for hexa h5
    update_new_hexahedron(d6, m_lcc->template create_attribute<3>(), nextlevel); // for hexa h6
    update_new_hexahedron(d7, m_lcc->template create_attribute<3>(), nextlevel); // for hexa h7

    if (split_volume_in_eight_functor()!=nullptr)
    { split_volume_in_eight_functor()(nfd); }
  }

public:
  // Get the dynamic splitedge in two functor.
  // This functor takes two darts in parameter, for the two edges resulting
  // of the split of the initial edge in two. This function is called after
  // the split, i.e. if 1-attributes are enabled, both edges have their own
  // 1-attribute correctly set. The function will only update the content of
  // the attributes.
  const std::function<void(Dart_descriptor, Dart_descriptor)>&
  split_edge_in_two_functor() const
  { return m_split_edge_in_two; }

  void set_split_edge_in_two_functor(std::function
                                     <void(Dart_descriptor, Dart_descriptor)> f)
  {
    m_split_edge_in_two_old = m_split_edge_in_two;
    m_split_edge_in_two = f;
  }

  void remove_split_edge_in_two_functor()
  { set_split_edge_in_two_functor(nullptr); }

  void restore_last_split_edge_in_two_functor()
  { set_split_edge_in_two_functor(m_split_edge_in_two_old); }

  // Get the dynamic splitface in four functor.
  // This functor takes four darts in parameter, for the four faces resulting
  // of the split of the initial face in four. Moreover, these dart belongs all
  // to one different edge incident to the middle vertex (and having this middle
  // vertex as extremity).
  // This function is called after the split, i.e. if 2-attributes are enabled,
  // the four faces have their own 2-attribute correctly set. The middle vertex
  // is also associated with a 0-attribute. However the 4 new edges are not (yet)
  // associated with 1-attributes.
  const std::function<void(Dart_descriptor, Dart_descriptor, Dart_descriptor, Dart_descriptor)>&
  split_face_in_four_functor() const
  { return m_split_face_in_four; }

  void set_split_face_in_four_functor(const std::function
                                     <void(Dart_descriptor, Dart_descriptor,
                                           Dart_descriptor, Dart_descriptor)>& f)
  {
    m_split_face_in_four_old = m_split_face_in_four;
    m_split_face_in_four = f;
  }

  void remove_split_face_in_four_functor()
  { set_split_face_in_four_functor(nullptr); }

  void restore_last_split_face_in_four_functor()
  { set_split_face_in_four_functor(m_split_face_in_four_old); }

  // Get the dynamic split volume in eight functor.
  // This functor takes an array of 12 dart handle as parameter:
  // First 4 darts: one for each right faces of the left volumes
  // Next 4 darts: one for each top faces of the low level volumes
  // Last 4 darts: one for each behind faces of the front volumes
  // Dart numbers 0, 1, 2, 3, 5, 7 allow to retreive the 6 new edges
  //          (incident to the central vertex)
  // Each dart in the array belongs to the central vertex
  const std::function<void(Dart_descriptor[])>&
  split_volume_in_eight_functor() const
  { return m_split_volume_in_eight; }

  void set_split_volume_in_eight_functor(const std::function
                                         <void(Dart_descriptor[])>& f)
  {
    m_split_volume_in_eight_old = m_split_volume_in_eight;
    m_split_volume_in_eight = f;
  }

  void remove_split_volume_in_eight_functor()
  { set_split_volume_in_eight_functor(nullptr); }

  void restore_last_split_volume_in_eight_functor()
  { set_split_volume_in_eight_functor(m_split_volume_in_eight_old); }


protected:
  /// The linear cell complex
  LCC* m_lcc;

  /// The mark for corner darts (dart that belong to a corner of a hexahedron)
  typename LCC::size_type m_corner_mark;

  /// The mark for edge direction: false = dir1; true = dir2
  /// @todo not generic, works only for square/hexahedra. For other type of
  /// elements/subdvision scheme, need a unsigned char for each dart have one type
  /// for each different type of edge => this will allow us to follow edge
  typename LCC::size_type m_edge_first_direction;

  // The subdvision level of each dart (that belongs to a hexahedron)
  CGAL::Unique_hash_map<Dart_descriptor, level_size,
  typename LCC::Hash_function> m_dart_levels;

  std::function<void(Dart_descriptor, Dart_descriptor)> m_split_edge_in_two;
  std::function<void(Dart_descriptor, Dart_descriptor)> m_split_edge_in_two_old;

  std::function<void(Dart_descriptor, Dart_descriptor, Dart_descriptor, Dart_descriptor)>
  m_split_face_in_four;

  std::function<void(Dart_descriptor, Dart_descriptor, Dart_descriptor, Dart_descriptor)>
  m_split_face_in_four_old;

  std::function<void(Dart_descriptor[])> m_split_volume_in_eight;
  std::function<void(Dart_descriptor[])> m_split_volume_in_eight_old;

  std::vector<Volume_info> m_volume_info;
};

#endif // HEXAHEDRAL_SUBDIVISION_H
