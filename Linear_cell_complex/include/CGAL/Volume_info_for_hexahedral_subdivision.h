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
#ifndef VOLUME_INFO_FOR_HEXAHEDRAL_SUBDIVISION_H
#define VOLUME_INFO_FOR_HEXAHEDRAL_SUBDIVISION_H

/**
 * Base class for an hexahedral 3D element in a LCC.
 * Each hexahedral stores its 6*4 corner darts, its three initial dimensions.
 *  Rule:
 *
 *       6----7
 *      /|   /|
 *     5----4 |
 *     | 3--|-2
 *     |/   |/
 *     0----1
 *
 * f0: 0, 5, 4, 1
 * f1: 1, 4, 7, 2
 * f2: 2, 7, 6, 3
 * f3: 3, 6, 5, 0
 * f4: 0, 1, 2, 3
 * f5: 4, 5, 6, 7
 *
 * e0:  05
 * e1:  14
 * e2:  27
 * e3:  36
 * e4:  01
 * e5:  12
 * e6:  23
 * e7:  30
 * e8:  45
 * e9:  56
 * e10: 67
 * e11: 74
 */
template<class Refs>
class Volume_info_for_hexahedral_subdivision
{
public:
  typedef typename Refs::Dart_descriptor DH;
  typedef typename Refs::FT FT;
  typedef typename Refs::Vertex_attribute_handle Vertex_attribute_handle;

  /*! Constructor */
  Volume_info_for_hexahedral_subdivision(unsigned char level=0) : m_subdivision_level(level)
  {
    // Dart corners
    for (unsigned int i=0; i<6; i++)
      for (unsigned int j=0; j<4; j++)
        m_corners[i][j] = NULL;
  }

  /*! Copy Constructor */
  Volume_info_for_hexahedral_subdivision(const Volume_info_for_hexahedral_subdivision& other) :
    m_subdivision_level(other.m_subdivision_level)
  {
    // Dart corners
    for (unsigned int i=0; i<6; i++)
      for (unsigned int j=0; j<4; j++)
       { m_corners[i][j] = other.m_corners[i][j]; }
  }

  /*! Operator = */
  void operator=(const Volume_info_for_hexahedral_subdivision& other)
  {
    m_subdivision_level = other.m_subdivision_level;
    // Dart corners
    for (unsigned int i=0; i<6; i++)
      for (unsigned int j=0; j<4; j++)
      { m_corners[i][j] = other.m_corners[i][j]; }
  }

  /*! To get corner dart -corner- of face -face- */
  DH get_corner_dart(unsigned int face, unsigned int corner)
  {
    assert(face<6 && corner<4);
    return m_corners[face][corner];
  }

  /*! To set corner dart -corner- of face -face- */
  void set_corner_dart(unsigned int face, unsigned int corner, DH dh)
  {
    assert(face<6 && corner<4);
    m_corners[face][corner]=dh;
  }

  /*! To get the -num_particle- th corner particle */
  template<typename LCC>
  Vertex_attribute_handle get_corner_particle(LCC& lcc, unsigned int num_particule)
  {
    assert(num_particule<8);

    // f4 : 0, 1, 2, 3
    if (num_particule<4)
    {
      return lcc.template attribute<0>(m_corners[4][num_particule]);
    }

    // f5 : 4, 5, 6, 7
    return lcc.template attribute<0>(m_corners[5][num_particule-4]);
  }

  /* ! To get a dart of the edge number e. */
  DH get_edge(unsigned int e)
  {
    assert(e<12);
    if (e<4)
    { // return the first corner dart in f[e]
      return m_corners[e][0];
    }
    else if (e<8)
    { // e=4,5,6 or 7: return the eth dart of face f4
      return m_corners[4][e-4];
    }
    // e=8,9,10 or 11 return the eth dart of face f5
    return m_corners[5][e-8];
  }

  /*! Initialize the corner dart of a specific face according to one dart */
  template<typename LCC>
  void initialize_face_corners (LCC& lcc, DH dh, unsigned int num_face)
  {
    assert(num_face<6);
    for (unsigned int i=0; i<4; i++)
    {
      m_corners[num_face][i] = dh;
      dh=lcc.template beta<1>(dh);
    }
  }

  /*! Initialization of the corner dart: 4 dart by face */
  template<typename LCC>
  void initialize_corners(LCC& lcc, DH dh)
  {
    // Dart corners
    for (unsigned int i=0; i<4; i++)
    {
      initialize_face_corners(lcc,dh,i);

      // Move to an other face
      dh=lcc.template beta<1,1,2>(dh);

    }//4 faces

    // Move to the down face
    dh=lcc.template beta<0,2>(dh);
    initialize_face_corners(lcc,dh,4);

    // Move to the up face
    dh=lcc.template beta<2,1,1,2>(dh);
    initialize_face_corners(lcc,dh,5);
  }

  DH (& get_corner_darts()) [6][4]
  { return m_corners; }

  unsigned char get_subdvision_level() const
  { return m_subdivision_level; }

  void set_subdivision_level(unsigned char l)
  { m_subdivision_level=l; }

  void increase_subdivision_level()
  {
    assert(m_subdivision_level<255);
    ++m_subdivision_level;
  }

  /*! To display information for debuging */
  template<typename LCC>
  void debug(LCC& lcc)
  {
    std::cout<<"Volume : "<<std::endl;
    for (unsigned int i=0; i<6; i++)
    {
      std::cout<<"  Face "<<i<<": ";
      for (unsigned int j=0; j<4; j++)
      {
        std::cout<<"dart: "<<&*(m_corners[i][j])<<" - "
                 <<lcc.point(m_corners[i][j])<<"; ";
      }
      std::cout<<std::endl;
    }
    std::cout<<"Corner particles: "<<std::endl;
    for (unsigned int i=0; i<8; i++)
    {
      std::cout<<"  Particle "<<i<<" : "
               <<&*(get_corner_particle(lcc,i))<<" - "
               <<lcc.point(get_corner_particle(lcc,i))<<std::endl;
    }
  }

protected:
  DH m_corners[6][4]; // 4 darts by face of the volume
  unsigned char m_subdivision_level;
};

#endif // VOLUME_INFO_FOR_HEXAHEDRAL_SUBDIVISION_H
