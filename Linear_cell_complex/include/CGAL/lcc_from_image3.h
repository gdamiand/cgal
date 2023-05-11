// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
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
#ifndef LCC_FROM_IMAGE3_H
#define LCC_FROM_IMAGE3_H
//******************************************************************************
#include <CGAL/Image_3.h>
#include <CGAL/Union_find.h>
#include <vector>
//******************************************************************************
namespace LCC_from_image_internal
{
template<typename LCC>
typename LCC::Dart_descriptor make_border(LCC& lcc, const CGAL::Image_3& im)
{
  typename LCC::Dart_descriptor first, last1, last2, tmp, tmp2;

  // 1. Create the first face.
  first=lcc.make_combinatorial_polygon(4);
  last1=lcc.beta(first, 1, 1);
  last2=lcc.beta(first, 1);

  // 2. We create the faces behind the first line of voxels.
  for(std::size_t x=0; x<=im.xdim(); ++x)
  {
    tmp=lcc.make_combinatorial_polygon(4);
    lcc.template link_beta<2>(tmp, last1);

    tmp2=lcc.make_combinatorial_polygon(4);
    lcc.template link_beta<2>(tmp2, last2);

    lcc.template link_beta<2>(lcc.beta(tmp, 1), lcc.beta(tmp2, 0));

    last1=lcc.beta(tmp, 1, 1);
    last2=lcc.beta(tmp2, 1, 1);
  }

  lcc.template link_beta<2>(first, last1);
  last1=lcc.beta(first, 1, 2, 1);

  // 3. Now we create the faces in the plane X-Y
  for(std::size_t y=1; y<=im.ydim(); ++y)
  {
    for(std::size_t x=0; x<=im.xdim(); ++x)
    {
      tmp=lcc.make_combinatorial_polygon(4);
      lcc.template link_beta<2>(lcc.beta(tmp, 0), last1);
      lcc.template link_beta<2>(tmp, last2);

      last1=lcc.beta(last1, 1, 2, 1);
      last2=lcc.beta(tmp, 1, 1);
    }
  }

  // 4. And now we link the borders.
  lcc.template link_beta<2>(last2, lcc.beta(first, 0));

  last2=lcc.beta(first, 1, 1, 2, 0);
  for(std::size_t x=0; x<=im.xdim(); ++x)
  {
    lcc.template link_beta<2>(last1, last2);
    last1=lcc.beta(last1, 1, 2, 1);
    last2=lcc.beta(last2, 0, 2, 0);
  }

  return lcc.beta(first, 1, 1);
}
//******************************************************************************
template<typename LCC>
void destroy_border(LCC& lcc, typename LCC::Dart_descriptor ADart,
                    const CGAL::Image_3& im)
{
  std::vector<typename LCC::Dart_descriptor> todelete;
  todelete.reserve(4*(2+im.xdim())*(2+im.ydim()));

  auto treated=lcc.get_new_mark();
  for(auto it=lcc.template darts_of_orbit_basic<1,2>(ADart, treated).begin(),
       itend=lcc.template darts_of_orbit_basic<1,2>(ADart, treated).end();
       it!=itend; ++it)
  { todelete.push_back(it); }

  for(auto dh: todelete)
  { lcc.erase_dart(dh); }

  assert(lcc.is_whole_map_unmarked(treated));
  lcc.free_mark(treated);
}
//******************************************************************************
template<typename LCC>
typename LCC::Dart_descriptor precode_l8(LCC& lcc,
                                         typename LCC::Dart_descriptor ALast,
                                         typename LCC::Dart_descriptor AUp,
                                         typename LCC::Dart_descriptor ABehind)
{
  assert(ALast!=NULL && AUp!=NULL && ABehind!=NULL);
  assert(!lcc.is_free(ALast, 2) && !lcc.is_free(AUp, 2) && !lcc.is_free(ABehind, 2));
  assert(lcc.beta(ALast, 2)==lcc.beta(ABehind, 0, 0));
  assert(lcc.beta(ABehind, 0, 2)==lcc.beta(AUp, 0, 0));
  assert(lcc.beta(ALast, 0, 2)==lcc.beta(AUp, 0));

  assert(lcc.beta(ALast, 1, 1, 1, 1)  ==ALast);
  assert(lcc.beta(AUp, 1, 1, 1, 1)    ==AUp);
  assert(lcc.beta(ABehind, 1, 1, 1, 1)==ABehind);

  typename LCC::Dart_descriptor t1=lcc.template beta<2>(AUp);
  typename LCC::Dart_descriptor t2=lcc.template beta<1, 2>(AUp);
  typename LCC::Dart_descriptor t3=lcc.template beta<1, 2>(ALast);
  typename LCC::Dart_descriptor t4=lcc.template beta<1, 2>(ABehind);

  // 1) First we modify the topology.

  // Left face => front face
  lcc.template link_beta<2>(lcc.beta(ALast, 0), t1);
  lcc.template link_beta<2>(lcc.beta(ALast, 1), AUp);

  // Behind face => right face
  lcc.template link_beta<2>(lcc.beta(ABehind, 0), t2);
  lcc.template link_beta<2>(lcc.beta(ABehind, 1), lcc.beta(AUp, 1));

  // Upper face => down face
  lcc.template link_beta<2>(lcc.beta(AUp, 0), t3);
  lcc.template link_beta<2>(lcc.beta(AUp, 1, 1), t4);

  return ABehind;
}
///////////////////////////////////////////////////////////////////////////////
template<typename Vector>
bool coplanar(const Vector& normal1, const Vector& normal2, double epsilon=0.1)
{
  double angle=CGAL::approximate_angle(normal1, normal2);
  return (angle<epsilon || angle>(360-epsilon) ||
          (angle>180-epsilon && angle<180+epsilon));
}
///////////////////////////////////////////////////////////////////////////////
/// Merge adjacent marked coplanar faces. Two faces are considered coplanar
/// if the angle between the two normal of the faces is smaller than epsilon.
/// A face is considered marked if one of its darts is marked.
template<class LCC>
void merge_coplanar_faces(LCC& lcc, typename LCC::size_type amark,
                          double epsilon=1.)
{
  lcc.set_update_attributes(false);

  typename LCC::size_type facemark=lcc.get_new_mark();

  /// Associative array that gives for each dart the index of its face
  std::unordered_map<typename LCC::Dart_handle, std::size_t> face_id;

  /// UF of faces to keep face connected
  typedef CGAL::Union_find<std::size_t> Face_UF;
  typedef Face_UF::handle Face_UF_handle;
  Face_UF faces_merged;

  /// Array of faces
  std::vector<std::pair<typename LCC::Vector, Face_UF_handle>> face_normal;

  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if(lcc.is_marked(it, amark) && !lcc.is_marked(it, facemark))
    {
      std::size_t nb=face_normal.size();
      face_normal.push_back
          (std::make_pair(CGAL::compute_normal_of_cell_2(lcc, it),
                          faces_merged.make_set(nb)));
      for(auto itd=lcc.template darts_of_cell_basic<2>(it, facemark).begin(),
          itdend=lcc.template darts_of_cell_basic<2>(it, facemark).end();
          itd!=itdend; ++itd)
      {
        lcc.mark(itd, facemark);
        face_id[itd]=nb;
      }
    }
  }

  std::vector<typename LCC::Dart_handle> edges_to_remove;
  for(auto f: face_id)
  {
    typename LCC::Dart_handle cur=f.first;
    do
    {
      if (lcc.is_marked(cur, facemark) &&
          lcc.is_marked(lcc.template beta<2>(cur), facemark) &&
          lcc.template is_removable<1>(cur) &&
          faces_merged.find(face_normal[f.second].second)!=
          faces_merged.find(face_normal[face_id[lcc.template beta<2>(cur)]].second))
      {
        if(coplanar(face_normal[f.second].first,
                    face_normal[face_id[lcc.template beta<2>(cur)]].first,
                    epsilon))
        {
          edges_to_remove.push_back(cur);
          faces_merged.unify_sets(face_normal[f.second].second,
              face_normal[face_id[lcc.template beta<2>(cur)]].second);
        }
      }
      cur=lcc.next(cur);
    }
    while(cur!=f.first);
  }

  for(auto e: edges_to_remove)
  {
    if (lcc.darts().is_used(e))
    { lcc.template remove_cell<1>(e); }
  }

  lcc.free_mark(facemark);
  lcc.set_update_attributes(true);
}
///////////////////////////////////////////////////////////////////////////////
/// Remove all marked dangling edges (and the edge becoming dangling after the
/// removal of dangling edges).
/// An edge is considered marked if one of its darts is marked.
template<class LCC>
void remove_danglig_edges(LCC& lcc, typename LCC::size_type amark)
{
  lcc.set_update_attributes(false);

  typename LCC::size_type edgemark=lcc.get_new_mark();
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if(lcc.is_marked(it, amark) && !lcc.is_marked(it, edgemark))
    { lcc.template mark_cell<1>(it, edgemark); }
  }

  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    typename LCC::Dart_handle cur=it, next;
    while(lcc.is_marked(cur, edgemark) && lcc.template beta<0, 2>(cur)==cur)
    {
      next=lcc.template beta<1>(cur);
      lcc.template remove_cell<1>(cur);
      cur=next;
    }
  }

  lcc.free_mark(edgemark);
  lcc.set_update_attributes(true);
}
//******************************************************************************
template<typename LCC>
void process_all_edges(LCC& lcc)
{
  auto amark=lcc.get_new_mark();
  lcc.negate_mark(amark); // mark all darts
  merge_coplanar_faces(lcc, amark);
  remove_danglig_edges(lcc, amark);
  lcc.free_mark(amark);
}
//******************************************************************************
template<typename LCC>
void process_all_vertices(LCC& lcc)
{
  typename LCC::size_type vertexmark=lcc.get_new_mark();
  lcc.negate_mark(vertexmark);
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if(lcc.is_marked(it, vertexmark))
    {
      if(lcc.template is_removable<0>(it)) { lcc.template remove_cell<0>(it); }
      else { lcc.template unmark_cell<0>(it, vertexmark); }
    }
  }

  lcc.free_mark(vertexmark);
}
//******************************************************************************
template<typename LCC>
typename LCC::Dart_descriptor create_cube(LCC& lcc,
                                          typename LCC::Dart_descriptor ALast,
                                          typename LCC::Dart_descriptor AUp,
                                          typename LCC::Dart_descriptor ABehind,
                                          const typename LCC::Point& p)
{
  typename LCC::Dart_descriptor res=lcc.make_combinatorial_hexahedron();

  // We need to use sew<3> due to the 3 additional vertex attribute creation.
  // Without, we could use topo_sew<3>. TODO check?
  lcc.template sew<3>(ALast, lcc.template beta<2, 1, 1>(res));
  lcc.template sew<3>(AUp, lcc.template beta<1, 2>(res));
  lcc.template sew<3>(ABehind, lcc.template beta<1, 1, 2, 1, 1, 2>(res));

  lcc.set_vertex_attribute(ALast, lcc.create_vertex_attribute(p));

  // The 3 following vertex attribute creations are not optimal, and maybe
  // useless, but maybe required in some disconnection cases... TODO check?
  lcc.set_vertex_attribute(lcc.template beta<2>(ALast),
                           lcc.create_vertex_attribute(typename LCC::Point(p.x(), p.y(), p.z()+1)));
  lcc.set_vertex_attribute(lcc.template beta<0, 2>(ABehind),
                           lcc.create_vertex_attribute(typename LCC::Point(p.x()+1, p.y(), p.z())));
  lcc.set_vertex_attribute(lcc.template beta<0, 2>(AUp),
                           lcc.create_vertex_attribute(typename LCC::Point(p.x(), p.y()+1, p.z())));

  return lcc.template beta<1, 1, 2, 1, 1>(res);
}
//******************************************************************************
template<typename LCC>
typename LCC::Dart_descriptor compute_up_from_last(LCC& lcc,
                                                   typename LCC::Dart_descriptor ADart)
{
  typename LCC::Dart_descriptor res=lcc.template beta<0, 2>(ADart);
  while(!lcc.template is_free<3>(res))
  { res=lcc.template beta<3, 2>(res); }
  return lcc.template beta<1>(res);
}
//******************************************************************************
template<typename LCC>
typename LCC::Dart_descriptor compute_behind_from_last(LCC& lcc,
                                                       typename LCC::Dart_descriptor ADart)
{
  typename LCC::Dart_descriptor res=lcc.template beta<2>(ADart);
  while(!lcc.template is_free<3>(res))
  { res=lcc.template beta<3, 2>(res); }
  return lcc.template beta<1, 1>(res);
}
}
//******************************************************************************
namespace CGAL
{
template<typename LCC>
bool lcc_from_image3(LCC& lcc, const CGAL::Image_3& im,
                     bool simplify_edges=false, bool simplify_vertices=false)
{
  using DD=typename LCC::Dart_descriptor;

  bool sameLeft;
  bool sameUp;
  bool sameBehind;

  DD up      =NULL;
  DD behind  =NULL;
  DD nextLast=NULL;

  typename LCC::Point p;

  DD last=LCC_from_image_internal::make_border(lcc, im);
  for(std::size_t z=0; z<=im.zdim(); ++z)
  {
    for(std::size_t y=0; y<=im.ydim(); ++y)
    {
      for(std::size_t x=0; x<=im.xdim(); ++x)
      {
        // 1. We compute the darts up and behind from the last dart.
        up    =LCC_from_image_internal::compute_up_from_last    (lcc, last);
        behind=LCC_from_image_internal::compute_behind_from_last(lcc, last);

        sameLeft=(y==im.ydim() || z==im.zdim()) ||
                   (x>0 && x<im.xdim() &&
                    im.value(x-1, y, z)==im.value(x, y, z));
        sameUp=(x==im.xdim() || y==im.ydim()) ||
                 (z>0 && z<im.zdim() &&
                  im.value(x, y, z-1)==im.value(x, y, z));
        sameBehind=(x==im.xdim() || z==im.zdim()) ||
                     (y>0 && y<im.ydim() &&
                      im.value(x, y-1, z)==im.value(x, y, z));

        p=typename LCC::Point(x, y, z);

        if(sameLeft && sameUp && sameBehind &&
            lcc.beta(last, 2)==lcc.beta(behind, 0, 0) &&
            lcc.beta(behind, 0, 2)==lcc.beta(up, 0, 0) &&
            lcc.beta(last, 0, 2)==lcc.beta(up, 0))
        {
          // Precode L8 is process independently for speed up.
          nextLast=LCC_from_image_internal::precode_l8(lcc, last, up, behind);
        }
        else
        {
          nextLast=LCC_from_image_internal::create_cube(lcc, last, up, behind, p);

          // We test the three possible face removal.
          if(sameLeft)
          { lcc.template remove_cell<2>(last); }

          if(sameBehind)
          { lcc.template remove_cell<2>(behind); }

          if(sameUp)
          { lcc.template remove_cell<2>(up); }

          /* std::cout<<"Voxel ("<<x<<", "<<y<<", "<<y<<"): ";
          lcc.display_characteristics(std::cout)<<" #points="<<lcc.number_of_vertex_attributes()<<std::endl; */
        }
        last = nextLast;
      }
    }
  }

  LCC_from_image_internal::destroy_border(lcc, last, im);

  if(simplify_edges)
  { LCC_from_image_internal::process_all_edges(lcc); }

  if(simplify_vertices)
  { LCC_from_image_internal::process_all_vertices(lcc); }

  lcc.display_characteristics(std::cout)<<" #points="<<lcc.number_of_vertex_attributes()<<std::endl;

  assert(lcc.is_valid());
  return true;
}
} // end namespace CGAL
//******************************************************************************
#endif // LCC_FROM_IMAGE3_H
